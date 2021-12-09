import sys
import math

from qgis.core import QgsGeometry, QgsFeature, QgsPoint, QgsPointXY,\
    QgsSpatialIndex, QgsFeatureRequest, QgsMessageLog, QgsGeometryUtils,\
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
    QgsMessageLog.logMessage(msg, "WaterQuickEdit")


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


def find_nearest_edge_from_multi_polygon(point_geometry, multipolygon_geometry):
    """
    @param point_geometry: QgsGeometry
    @param multipolygon_geometry: QgsGeometry
    @return list of QgsPointXY
    """

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
            log_message(str(distance))

            min_distance = distance
            min_start_point = start_point
            min_end_point = end_point

    return [min_start_point, min_end_point]


def find_nearest_line_segment(point_geometry,  multipolyline_feature):
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
    length = (start_p.x() - end_p.x()) * (start_p.x() - end_p.x()) + (start_p.y() - end_p.y()) * (start_p.y() - end_p.y())

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







