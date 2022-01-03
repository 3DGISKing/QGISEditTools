# -*- coding: utf-8 -*-

import os
import math

# noinspection PyUnresolvedReferences
from PyQt5.QtCore import pyqtSlot, QSettings, QDate, QLocale, QTranslator, QCoreApplication, Qt, QTimer, QVariant

# noinspection PyUnresolvedReferences
from PyQt5.QtGui import (QIcon, QColor, QCursor, QPixmap)

# noinspection PyUnresolvedReferences
from PyQt5.QtWidgets import (QAction, QApplication, QDialog, QInputDialog, QStyle, QMenu, QMenuBar, QMessageBox, QToolBar,
                             QDockWidget)

# noinspection PyUnresolvedReferences
from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsMessageLog, Qgis, QgsPointXY, QgsProject, \
    QgsWkbTypes, QgsPoint, QgsVectorLayerEditUtils

# noinspection PyUnresolvedReferences
from qgis.gui import QgsHighlight, QgsRubberBand, QgsMapToolEmitPoint, QgsVertexMarker

# noinspection PyUnresolvedReferences
from PyQt5 import QtGui, QtWidgets

from qgis.core import QgsWkbTypes
from qgis.core import QgsVectorLayer
from qgis.core import QgsProject, QgsMessageLog, QgsFeature, QgsGeometry, QgsField, QgsVectorLayerUtils,QgsRectangle
from .util import get_cut_distance_by_area_limit, cut_polygon, distance_between_two_points, find_nearest_feature_from_features, find_nearest_edge_from_multi_polygon,find_perpendicular_point, generate_intersects, log_message


class QGISEditTools:
    def __init__(self, iface):
        # Save reference to the QGIS interface
        self._iface = iface
        self.map_canvas = iface.mapCanvas()
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        self._divide_by_distance_button_clicked = False
        self._divide_by_area_button_clicked = False

    # noinspection PyAttributeOutsideInit
    # noinspection PyPep8Naming
    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, "icons", "divide.png")
        self._action_polygon_divide_by_distance = QAction(QIcon(icon_path), 'Divide by distance', self._iface.mainWindow())
        self._action_polygon_divide_by_area = QAction(QIcon(icon_path), 'Divide by area', self._iface.mainWindow())

        self._toolbar = self._iface.addToolBar(u'QGIS Edit')

        self._toolbar.addAction(self._action_polygon_divide_by_distance)
        self._toolbar.addAction(self._action_polygon_divide_by_area)

        self.reference_line_selector = QgsMapToolEmitPoint(self.map_canvas)

        reference_line_marker = QgsVertexMarker(self.map_canvas)

        reference_line_marker.setIconSize(10)
        reference_line_marker.setIconType(QgsVertexMarker.ICON_BOX)
        reference_line_marker.setPenWidth(3)
        reference_line_marker.setColor(QColor(0, 255, 0))

        self.reference_line_marker = reference_line_marker

        line_rubber_band = QgsRubberBand(self.map_canvas, QgsWkbTypes.LineGeometry)
        line_rubber_band.setStrokeColor(QColor('red'))
        line_rubber_band.setWidth(3)

        self.line_rubber_band = line_rubber_band

        line_rubber_band = QgsRubberBand(self.map_canvas, QgsWkbTypes.LineGeometry)
        line_rubber_band.setStrokeColor(QColor('red'))
        line_rubber_band.setWidth(3)

        self.cut_line_rubber_band = line_rubber_band

        self._show_cut_line = False

        # connecting
        self._action_polygon_divide_by_distance.triggered.connect(self.on_action_divide_by_distance_triggered)
        self.reference_line_selector.canvasClicked.connect(self.on_reference_line_selector_canvas_clicked)

        self._action_polygon_divide_by_area.triggered.connect(self.on_action_divide_by_area)

        self._flash_timer = QTimer(self.map_canvas)
        self._flash_timer.timeout.connect(self._on_flash_timer_timeout)

        self._highlight_list = []

    def unload(self):
        self.map_canvas.scene().removeItem(self.reference_line_marker)
        self.map_canvas.scene().removeItem(self.line_rubber_band)
        self.map_canvas.scene().removeItem(self.cut_line_rubber_band)

        del self._toolbar

    def check_polygon_divide_condition(self):
        layers = self._iface.layerTreeView().selectedLayers()

        if len(layers) == 0:
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'No selected layer!')
            return False

        active_layer = self._iface.activeLayer()

        if active_layer is None:
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'No active layer!')
            return False

        if active_layer.__class__.__name__ != 'QgsVectorLayer':
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'Selected layer is not Vector layer!')
            return False

        geom_type = active_layer.geometryType()

        if geom_type is not None:
            if geom_type != QgsWkbTypes.PolygonGeometry:
                QMessageBox.information(self._iface.mainWindow(), 'Info', 'Selected layer is not Polygon layer!')
                return False

        selected_features = active_layer.selectedFeatures()

        if len(selected_features) == 0:
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "No selected")
            return False

        return True

    def on_action_divide_by_distance_triggered(self):
        if not self.check_polygon_divide_condition():
            return

        self._divide_by_distance_button_clicked = True
        self.map_canvas.setMapTool(self.reference_line_selector)

    def on_action_reference_line_selector_triggered(self):
        self.map_canvas.setMapTool(self.reference_line_selector)

    def hide_all_marker_rubber_band(self):
        self.reference_line_marker.hide()
        self.line_rubber_band.hide()

    def mark_reference_line(self, edge):
        start_point = edge[0]
        end_point = edge[1]

        # get mid point of the nearest edge
        mid_x = (start_point.x() + end_point.x()) / 2
        mid_y = (start_point.y() + end_point.y()) / 2

        self.reference_line_marker.setCenter(QgsPointXY(mid_x, mid_y))
        self.reference_line_marker.show()

        polyline = QgsGeometry.fromPolylineXY([start_point, end_point])

        self.line_rubber_band.setToGeometry(polyline)
        self.line_rubber_band.show()

    def on_reference_line_selector_canvas_clicked(self, point, button):
        """
        :type: QgsPointXY
        :param: point

        """
        if self._divide_by_distance_button_clicked:
            self.divide_by_distance(point)
        elif self._divide_by_area_button_clicked:
            self.divide_by_area(point)

    def divide_by_distance(self, point):
        active_layer = self._iface.activeLayer()
        selected_features = active_layer.selectedFeatures()
        nearest_feature = find_nearest_feature_from_features(QgsGeometry.fromPointXY(point), selected_features)

        # get start and end of nearest edge of the nearest polygon
        nearest_edge = find_nearest_edge_from_multi_polygon(point, selected_features)

        self.mark_reference_line(nearest_edge)
        distance, ok = QInputDialog.getDouble(self._iface.mainWindow(), "Input", "enter a distance")

        if not ok:
            self.hide_all_marker_rubber_band()
            return

        self.do_cut_polygon(nearest_feature, nearest_edge, distance)
        self._divide_by_distance_button_clicked = False

    def divide_by_area(self, point):
        active_layer = self._iface.activeLayer()
        selected_features = active_layer.selectedFeatures()
        nearest_feature = find_nearest_feature_from_features(QgsGeometry.fromPointXY(point), selected_features)

        # get start and end of nearest edge of the nearest polygon
        nearest_edge = find_nearest_edge_from_multi_polygon(point, selected_features)

        self.mark_reference_line(nearest_edge)

        area, ok = QInputDialog.getDouble(self._iface.mainWindow(), "Input", "enter a area")

        if area < 0:
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "Invalid area!")
            return

        if area > nearest_feature.geometry().area():
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "Invalid area!. Selected feature 's area is " +
                                 str(nearest_feature.geometry().area()))
            return

        if not ok:
            self.hide_all_marker_rubber_band()
            return

        distance = get_cut_distance_by_area_limit(nearest_feature, nearest_edge[0], nearest_edge[1], area)
        self.do_cut_polygon(nearest_feature, nearest_edge, distance)
        self._divide_by_area_button_clicked = False

    def do_cut_polygon(self, nearest_feature, nearest_edge, distance):
        geometry = nearest_feature.geometry()

        success, split_geometry_list, topo, cut_line = cut_polygon(geometry, nearest_edge[0],
                                                                   nearest_edge[1],
                                                                   distance)

        if self._show_cut_line:
            self.cut_line_rubber_band.setToGeometry(cut_line)
            self.cut_line_rubber_band.show()

        if success != QgsGeometry.OperationResult.Success:
            self.hide_all_marker_rubber_band()
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "failed to split")
            return

        if len(split_geometry_list) != 1:
            self.hide_all_marker_rubber_band()
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "failed to split")
            return

        active_layer = self._iface.activeLayer()
        active_layer.startEditing()

        # selected feature's geometry changed
        active_layer.changeGeometry(nearest_feature.id(), geometry)

        # add new feature
        feature = QgsVectorLayerUtils.createFeature(active_layer)

        feature.setGeometry(split_geometry_list[0])
        active_layer.addFeature(feature)

        active_layer.commitChanges()

        active_layer.removeSelection()

        # feature_ids = [nearest_feature.id(), feature.id()]
        # self.map_canvas.flashFeatureIds(active_layer, feature_ids)

        self.hide_all_marker_rubber_band()
        self._flash_geometries(active_layer, [nearest_feature.geometry(), split_geometry_list[0]])
        # self._flash_geometries(active_layer, [split_geometry_list[0]])

    def _on_flash_timer_timeout(self):
        self._flash_timer.stop()

        for item in self._highlight_list:
            self.map_canvas.scene().removeItem(item)

        self._highlight_list = []

    def _flash_geometries(self, layer, geometries):
        for geometry in geometries:
            h = QgsHighlight(self.map_canvas, geometry, layer)
            h.setColor(QColor(255, 0, 0, 255))
            h.setWidth(5)
            h.setFillColor(QColor(255, 0, 0, 100))
            self._highlight_list.append(h)

        self._flash_timer.start(2000)  # Milliseconds before finishing the flash

    def on_action_divide_by_area(self):
        if not self.check_polygon_divide_condition():
            return

        self._divide_by_area_button_clicked = True
        self.map_canvas.setMapTool(self.reference_line_selector)


if __name__ == '__main__':
    pass
