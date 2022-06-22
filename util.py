import sys
import math

from qgis.core import QgsGeometry, QgsFeature, QgsPoint, QgsPointXY, \
    QgsSpatialIndex, QgsFeatureRequest, QgsMessageLog, QgsGeometryUtils, \
    QgsVectorLayer, QgsWkbTypes

RADIUS = 0.25


def calc_intersection(
        start_point_line_a,  # type: QgsPoint
        end_point_line_a,  # type: QgsPoint
        start_point_line_b,  # type: QgsPoint
        end_point_line_b  # type: QgsPoint
):
    line_a = QgsGeometry.fromPolylineXY([start_point_line_a, end_point_line_a])
    line_b = QgsGeometry.fromPolylineXY([start_point_line_b, end_point_line_b])

    intersection_geometry = line_a.intersection(line_b)
    intersection_point = intersection_geometry.asPoint()

    return intersection_point


def create_polygon(point_array, first_attribute_value):
    polygon_geometry = QgsGeometry.fromPolygonXY([point_array])

    feature = QgsFeature()
    feature.setGeometry(polygon_geometry)

    # prepare one attribute
    feature.initAttributes(1)
    feature.setAttribute(0, first_attribute_value)

    return feature


def divide(polygon_feature, num_1000_value):
    # polygon = polygon_feature.geometry().asPolygon()
    multi_polygon = polygon_feature.geometry().asMultiPolygon()

    polygon = multi_polygon[0]

    n = len(polygon[0])
    mid_point_array = []

    for i in range(n):
        point = polygon[0][i]

        if i < 3:
            next_point = polygon[0][i + 1]
        else:
            next_point = polygon[0][0]

        x = point.x()
        y = point.y()

        next_x = next_point.x()
        next_y = next_point.y()

        mid_x = (x + next_x) / 2
        mid_y = (y + next_y) / 2

        mid_point_array.append(QgsPointXY(mid_x, mid_y))

        # QgsMessageLog.logMessage(str(x) + ' ' + str(y), 'IROOMGIS_Tools')

    center = calc_intersection(mid_point_array[0], mid_point_array[2], mid_point_array[1], mid_point_array[3])

    # QgsMessageLog.logMessage(num_10_value, 'IROOMGIS_Tools')

    divided_polygons = []

    points = [polygon[0][0], mid_point_array[0], center, mid_point_array[3]]
    divided_polygons.append(create_polygon(points, num_1000_value + "A"))

    points = [mid_point_array[0], polygon[0][1], mid_point_array[1], center]
    divided_polygons.append(create_polygon(points, num_1000_value + "B"))

    points = [center, mid_point_array[1], polygon[0][2], mid_point_array[2]]
    divided_polygons.append(create_polygon(points, num_1000_value + "D"))

    points = [mid_point_array[3], center, mid_point_array[2], polygon[0][3]]
    divided_polygons.append(create_polygon(points, num_1000_value + "C"))

    return divided_polygons


def log_message(msg):
    QgsMessageLog.logMessage(msg, "QGISEditTool")


def find_perpendicular_point(point_feature, line_layer):
    point_geometry = point_feature.geometry()
    nearest_line_feature = find_nearest_feature(point_geometry, line_layer)
    projected_point = find_nearest_line_segment(point_geometry, nearest_line_feature)

    return projected_point


def find_nearest_feature_old(
        point_geometry,  # type: QgsPointXY
        layer  # type: QgsVectorLayer
):
    provider = layer.dataProvider()
    spatial_index = QgsSpatialIndex()  # create spatial index object

    feature = QgsFeature()
    features = provider.getFeatures()  # gets all features in layer

    # insert features to index
    while features.nextFeature(feature):
        spatial_index.insertFeature(feature)

    nearest_ids = spatial_index.nearestNeighbor(point_geometry, 1)  # we need only one neighbour

    feature_id = nearest_ids[0]

    ids = [feature_id]

    features = layer.getFeatures(QgsFeatureRequest().setFilterFid(feature_id))

    features.nextFeature(feature)

    return feature


def find_nearest_feature_from_features(geometry, features):
    """
    :param geometry:
    :type geometry: QgsGeometry

    :param features:
    :type features: list

    :return:
    :rtype: QgsFeature
    """

    min_distance = sys.float_info.max
    nearest_feature = None

    for feature in features:
        dist = feature.geometry().distance(geometry)
        if dist < min_distance:
            min_distance = dist
            nearest_feature = feature
    return nearest_feature


def find_nearest_feature(geometry, layer):
    """
    :param geometry:
    :type geometry: QgsGeometry

    :param layer:
    :type layer: QgsVectorLayer

    :return:
    :rtype: QgsFeature
    """
    provider = layer.dataProvider()
    features = provider.getFeatures()  # gets all features in layer

    min_distance = sys.float_info.max
    nearest_feature = None

    for feature in features:
        dist = feature.geometry().distance(geometry)

        if dist < min_distance:
            min_distance = dist
            nearest_feature = feature

    return nearest_feature


def find_nearest_edge_from_multi_polygon(point, features):
    """
    @param point: QgsGeometry
    @param features: QgsGeometry
    @return list of QgsPointXY
    """

    point_geometry = QgsGeometry.fromPointXY(point)
    nearest_feature = find_nearest_feature_from_features(point_geometry, features)
    multipolygon_geometry = nearest_feature.geometry()

    assert multipolygon_geometry.wkbType() == QgsWkbTypes.MultiPolygon, "error"

    multipolygon = multipolygon_geometry.asMultiPolygon()

    first_multipolygon = multipolygon[0]

    exterior_ring = first_multipolygon[0]

    min_distance = sys.float_info.max

    min_start_point = None
    min_end_point = None

    for i in range(len(exterior_ring) - 1):
        start_point = exterior_ring[i]
        end_point = exterior_ring[i + 1]

        points = [start_point, end_point]

        line_segment_geometry = QgsGeometry.fromPolylineXY(points)

        distance = line_segment_geometry.distance(point_geometry)

        if distance < min_distance:
            min_distance = distance
            min_start_point = start_point
            min_end_point = end_point

    return [min_start_point, min_end_point]


def find_nearest_line_segment(point_geometry, multipolyline_feature):
    """
    :param point_geometry:
    :type point_geometry: QgsPointXY

    :param multipolyline_feature:
    :type multipolyline_feature QgsFeature

    :return:
    :rtype: QgsPoint
    """
    geometry = multipolyline_feature.geometry()
    multi_polyline = geometry.asMultiPolyline()
    first_polyline = multi_polyline[0]

    min_distance = sys.float_info.max
    min_start_point = None
    min_end_point = None

    for i in range(len(first_polyline) - 1):
        start_point = first_polyline[i]
        end_point = first_polyline[i + 1]

        points = [start_point, end_point]

        # assume point is QgsPointXY
        line_segment = QgsGeometry.fromPolylineXY(points)

        distance = line_segment.distance(point_geometry)

        if distance < min_distance:
            min_distance = distance
            min_start_point = start_point
            min_end_point = end_point

    p = point_geometry.asPoint()
    p = QgsPoint(p)  # p is QgsPointXY

    s1 = QgsPoint(min_start_point)
    s2 = QgsPoint(min_end_point)

    projected_point = QgsGeometryUtils.projectPointOnSegment(p, s1, s2)

    return projected_point


def generate_intersects(feature1, feature2):
    geometry1 = feature1.geometry()
    geometry2 = feature2.geometry()

    # need to check wkb type
    # for now we assume multi polyline
    # geometry1.wkbType()

    # log_message("wkbType: " + str(geometry1.wkbType()))

    multi_polyline1 = geometry1.asMultiPolyline()
    multi_polyline2 = geometry2.asMultiPolyline()

    for i in range(len(multi_polyline1)):
        polyline1 = multi_polyline1[i]
        polyline1_geom = QgsGeometry.fromPolylineXY(polyline1)

        for j in range(len(multi_polyline2)):
            polyline2 = multi_polyline2[j]

            polyline2_geom = QgsGeometry.fromPolylineXY(polyline2)

            if polyline1_geom.intersects(polyline2_geom):
                result_polyline = generate_intersects_points(polyline1, polyline2)
                return result_polyline

    intersection_geom = feature1.geometry().intersection(feature2.geometry())
    intersection_point = intersection_geom.asPoint()

    log_message("2 " + str(intersection_point.x()) + " " + str(intersection_point.y()))


def distance_between_two_points(start_p, end_p):
    length = (start_p.x() - end_p.x()) * (start_p.x() - end_p.x()) + (start_p.y() - end_p.y()) * (
                start_p.y() - end_p.y())

    return math.sqrt(length)


def calc_rotated_position(p_x, p_y, pivot_x, pivot_y, angle):
    x = p_x - pivot_x
    y = p_y - pivot_y

    xx = math.cos(angle) * x - math.sin(angle) * y
    yy = math.sin(angle) * x + math.cos(angle) * y

    x = xx + pivot_x
    y = yy + pivot_y

    return [x, y]


def generate_intersects_points(polyline_xy1, polyline_xy2):
    result_polyline = []

    for i in range(len(polyline_xy1) - 1):
        start_point_1 = polyline_xy1[i]
        end_point_1 = polyline_xy1[i + 1]

        line_segment_geometry1 = QgsGeometry.fromPolylineXY([start_point_1, end_point_1])

        for j in range(len(polyline_xy2) - 1):
            start_point_2 = polyline_xy2[j]
            end_point_2 = polyline_xy2[j + 1]

            line_segment_geometry2 = QgsGeometry.fromPolylineXY([start_point_2, end_point_2])

            if not line_segment_geometry1.intersects(line_segment_geometry2):
                continue

            distance = distance_between_two_points(start_point_2, end_point_2)

            if distance < RADIUS * 2:
                continue

            intersection_geom = line_segment_geometry1.intersection(line_segment_geometry2)
            intersection_point = intersection_geom.asPoint()
            distance_from_start_to_intersection_2 = distance_between_two_points(start_point_2, intersection_point)

            radius_start_p_x = start_point_2.x() + (distance_from_start_to_intersection_2 - RADIUS) * (
                    end_point_2.x() - start_point_2.x()) / distance
            radius_start_p_y = start_point_2.y() + (distance_from_start_to_intersection_2 - RADIUS) * (
                    end_point_2.y() - start_point_2.y()) / distance

            radius_end_p_x = start_point_2.x() + (distance_from_start_to_intersection_2 + RADIUS) * (
                    end_point_2.x() - start_point_2.x()) / distance
            radius_end_p_y = start_point_2.y() + (distance_from_start_to_intersection_2 + RADIUS) * (
                    end_point_2.y() - start_point_2.y()) / distance

            for k in range(j + 1):
                point = polyline_xy2[k]
                result_polyline.append(point)

            # result_polyline.append(start_point_2)
            result_polyline.append(QgsPointXY(radius_start_p_x, radius_start_p_y))

            angle = math.pi / 8
            sign = 1

            rotated_position = calc_rotated_position(radius_start_p_x, radius_start_p_y, intersection_point.x(),
                                                     intersection_point.y(), angle)
            # compare y
            if rotated_position[1] < radius_start_p_y:
                sign = -1

            for k in range(1, 8):
                angle = math.pi / 8 * k * sign

                rotated_position = calc_rotated_position(radius_start_p_x, radius_start_p_y, intersection_point.x(),
                                                         intersection_point.y(), angle)

                result_polyline.append(QgsPointXY(rotated_position[0], rotated_position[1]))

            result_polyline.append(QgsPointXY(radius_end_p_x, radius_end_p_y))
            # result_polyline.append(end_point_2)

            for k in range(j + 1, len(polyline_xy2)):
                point = polyline_xy2[k]
                result_polyline.append(point)

            return result_polyline


def cut_polygon(polygon_geometry, start_point_of_edge, end_point_edge, distance):
    distance = -distance

    edge = QgsGeometry.fromPolylineXY([start_point_of_edge, end_point_edge])

    offset_edge = edge.offsetCurve(distance, 8, QgsGeometry.JoinStyleRound, 5.0)

    offset_edge_start = offset_edge.asPolyline()[0]
    offset_edge_end = offset_edge.asPolyline()[1]

    bounding_box = polygon_geometry.boundingBox()

    extend_length = max(bounding_box.width(), bounding_box.height())

    edge_length = distance_between_two_points(start_point_of_edge, end_point_edge)

    # get mid point of the offset_edge
    mid_x = (offset_edge_start.x() + offset_edge_end.x()) / 2
    mid_y = (offset_edge_start.y() + offset_edge_end.y()) / 2

    # extend offset_edge by edge_length / 2
    new_start_x = offset_edge_start.x() + (offset_edge_start.x() - mid_x) / (edge_length / 2) * extend_length
    new_start_y = offset_edge_start.y() + (offset_edge_start.y() - mid_y) / (edge_length / 2) * extend_length

    new_end_x = offset_edge_end.x() + (offset_edge_end.x() - mid_x) / (edge_length / 2) * extend_length
    new_end_y = offset_edge_end.y() + (offset_edge_end.y() - mid_y) / (edge_length / 2) * extend_length

    cut_line = QgsGeometry.fromPolylineXY([QgsPointXY(new_start_x, new_start_y), QgsPointXY(new_end_x, new_end_y)])

    #  success, split_geometry_list, topo
    ret = polygon_geometry.splitGeometry(cut_line.asPolyline(), False)

    # tuple add
    ret = ret + (cut_line,)

    return ret


# check divided polygon area is smaller than area limit

def cut_polygon_area_le_limit(polygon_geometry, start_point_of_edge, end_point_edge, distance, area_limit):
    distance = -distance

    success, split_geometry_list, topo, cut_line = cut_polygon(polygon_geometry, start_point_of_edge, end_point_edge,
                                                               distance)

    if success != QgsGeometry.OperationResult.Success:
        log_message("no success")
        return False

    if len(split_geometry_list) != 1:
        log_message("len invalid")
        return False

    area = split_geometry_list[0].area()

    log_message("area: " + str(area) + " at distance :" + str(-distance))

    return area <= area_limit


TOLERANCE = 0.01


def get_cut_distance_by_area_limit(polygon_feature, start_point_of_edge, end_point_edge, area_limit):
    wkt = polygon_feature.geometry().asWkt()
    geometry = QgsGeometry.fromWkt(wkt)
    bounding_box = geometry.boundingBox()
    max_area = bounding_box.width() * bounding_box.height()
    min_width = area_limit / bounding_box.height()

    distance_step = 0.1
    distance = distance_step
    area_error = area_limit

    while area_error > TOLERANCE:
        geometry = QgsGeometry.fromWkt(wkt)
        success, split_geometry_list, topo, cut_line = cut_polygon(geometry, start_point_of_edge, end_point_edge,
                                                                   distance)

        if success != QgsGeometry.OperationResult.Success:
            return -1

        if len(split_geometry_list) != 1:
            return -1

        area = split_geometry_list[0].area()

        while area <= area_limit:
            distance += distance_step
            geometry = QgsGeometry.fromWkt(wkt)
            success, split_geometry_list, topo, cut_line = cut_polygon(geometry, start_point_of_edge, end_point_edge,
                                                                       distance)

            if success != QgsGeometry.OperationResult.Success:
                return -1

            if len(split_geometry_list) != 1:
                return -1

            area = split_geometry_list[0].area()

        area_error = abs(area_limit - split_geometry_list[0].area())

        log_message("area: " + str(area) + " error: " + str(area_error) + " at distance : " + str(distance))

        # go back one step
        distance = distance - distance_step
        distance_step = distance_step / 10.0

    return distance

# remove vertice based on half plane divided by edge

def get_fixed_polygon(geometry, edge, ignore_minus_half_plane):
    assert geometry.wkbType() == QgsWkbTypes.Polygon

    start_point = edge[0]
    end_point = edge[1]

    x1 = start_point.x()
    y1 = start_point.y()

    x2 = end_point.x()
    y2 = end_point.y()

    # find the equation of line: Ax + By + C = 0

    a = (y2 - y1) / (x2 - x1)
    b = -1
    c = y1 - a * x1

    polygon = geometry.asPolygon()

    fixed_polygon = []
    for j in range(len(polygon)):
        ring = polygon[j]

        log_message(str(j) + " ring point, count:" + str(len(ring)) + "\n")

        fixed_ring = []

        for k in range(len(ring)):
            point = ring[k]
            x = point.x()
            y = point.y()

            # calculate the signed distance from point to edge
            signed_distance = a * x + b * y + c
            signed_distance = signed_distance / math.sqrt(a * a + b * b)

            if abs(signed_distance) < 0.1:
                log_message("x = " + str(x) + ', y = ' + str(y) + " aligned with edge")
                fixed_ring.append(point)
            else:
                if signed_distance > 0:
                    log_message("x = " + str(x) + ', y = ' + str(y) + "sign +")

                    if ignore_minus_half_plane:
                        fixed_ring.append(point)

                else:
                    if not ignore_minus_half_plane:
                        fixed_ring.append(point)

                    log_message("x = " + str(x) + ', y = ' + str(y) + "sign -")

        if len(fixed_ring) > 0:
            fixed_polygon.append(fixed_ring)

    return QgsGeometry.fromPolygonXY(fixed_polygon)

def get_fixed_multi_polygon(geometry, edge, ignore_minus_half_plane):
    assert geometry.wkbType() == QgsWkbTypes.MultiPolygon

    start_point = edge[0]
    end_point = edge[1]

    x1 = start_point.x()
    y1 = start_point.y()

    x2 = end_point.x()
    y2 = end_point.y()

    # find the equation of line: Ax + By + C = 0

    a = (y2 - y1) / (x2 - x1)
    b = -1
    c = y1 - a * x1

    fixed_multi_polygon = []

    multi_polygon = geometry.asMultiPolygon()

    polygon_count = len(multi_polygon)

    log_message("polygon count " + str(polygon_count))

    for i in range(polygon_count):
        polygon = multi_polygon[i]

        area = QgsGeometry.fromPolygonXY(polygon).area()

        if area < 0.1:
            log_message(str(i) + "th polygon area " + str(area) + " ignored")
            continue

        fixed_polygon = []
        for j in range(len(polygon)):
            ring = polygon[j]

            log_message(str(i) + "th polygon " + str(j) + " ring point, count:" + str(len(ring)) + "\n")

            fixed_ring = []

            for k in range(len(ring)):
                point = ring[k]
                x = point.x()
                y = point.y()

                # calculate the signed distance from point to edge
                signed_distance = a * x + b * y + c
                signed_distance = signed_distance / math.sqrt(a * a + b * b)

                if abs(signed_distance) < 0.1:
                    log_message("x = " + str(x) + ', y = ' + str(y) + " aligned with edge")
                    fixed_ring.append(point)
                else:
                    if signed_distance > 0:
                        log_message("x = " + str(x) + ', y = ' + str(y) + "sign +")
                        if ignore_minus_half_plane:
                            fixed_ring.append(point)
                    else:
                        log_message("x = " + str(x) + ', y = ' + str(y) + "sign -")
                        if not ignore_minus_half_plane:
                            fixed_ring.append(point)

            if len(fixed_ring) > 0:
                fixed_polygon.append(fixed_ring)

        if len(fixed_polygon) > 0:
            fixed_multi_polygon.append(fixed_polygon)

    return QgsGeometry.fromMultiPolygonXY(fixed_multi_polygon)

def get_fixed_polygon_geometry(geometry, edge, ignore_minus_half_plane):
    if geometry.wkbType() == QgsWkbTypes.MultiPolygon:
        return get_fixed_multi_polygon(geometry, edge, ignore_minus_half_plane)
    elif geometry.wkbType() == QgsWkbTypes.Polygon:
        return get_fixed_polygon(geometry, edge, ignore_minus_half_plane)
    elif geometry.wkbType() == QgsWkbTypes.GeometryCollection:
        geometry_collection = geometry.asGeometryCollection()

        for i in range(len(geometry_collection)):
            sub_geometry = geometry_collection[i]

            if sub_geometry.wkbType() == QgsWkbTypes.Polygon:
                fixed_geom = get_fixed_polygon(sub_geometry, edge, ignore_minus_half_plane)

                if fixed_geom is not None:
                    return fixed_geom
            elif sub_geometry.wkbType() == QgsWkbTypes.MultiPolygon:
                fixed_geom = get_fixed_multi_polygon(sub_geometry, edge, ignore_minus_half_plane)

                if fixed_geom is not None:
                    return fixed_geom