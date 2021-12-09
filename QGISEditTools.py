# -*- coding: utf-8 -*-

import os
import math

# noinspection PyUnresolvedReferences
from PyQt5.QtCore import pyqtSlot, QSettings, QDate, QLocale, QTranslator, QCoreApplication, Qt, QVariant

# noinspection PyUnresolvedReferences
from PyQt5.QtGui import (QIcon, QColor, QCursor, QPixmap)

# noinspection PyUnresolvedReferences
from PyQt5.QtWidgets import (QAction, QApplication, QDialog, QInputDialog, QStyle, QMenu, QMenuBar, QMessageBox, QToolBar,
                             QDockWidget)

# noinspection PyUnresolvedReferences
from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsMessageLog, Qgis, QgsPointXY, QgsProject, \
    QgsWkbTypes, QgsPoint, QgsVectorLayerEditUtils

# noinspection PyUnresolvedReferences
from qgis.gui import QgsRubberBand, QgsMapToolEmitPoint, QgsVertexMarker

# noinspection PyUnresolvedReferences
from PyQt5 import QtGui, QtWidgets

from qgis.core import QgsWkbTypes
from qgis.core import QgsVectorLayer
from qgis.core import QgsProject, QgsMessageLog, QgsFeature, QgsGeometry, QgsField, QgsVectorLayerUtils

from .util import distance_between_two_points, find_nearest_feature, find_nearest_edge_from_multi_polygon,find_perpendicular_point, generate_intersects, log_message


class QGISEditTools:
    def __init__(self, iface):
        # Save reference to the QGIS interface
        self._iface = iface
        self.map_canvas = iface.mapCanvas()
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

    # noinspection PyAttributeOutsideInit
    # noinspection PyPep8Naming
    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, "icons", "divide.png")
        self._action_divide = QAction(QIcon(icon_path), 'Divide', self._iface.mainWindow())

        self._toolbar = self._iface.addToolBar(u'IROOMGIS_Tools')

        self._toolbar.addAction(self._action_divide)

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
        self._action_divide.triggered.connect(self.on_action_divide_triggered)
        self.reference_line_selector.canvasClicked.connect(self.on_reference_line_selector_canvas_clicked)

    def unload(self):
        self.map_canvas.scene().removeItem(self.reference_line_marker)
        self.map_canvas.scene().removeItem(self.line_rubber_band)
        self.map_canvas.scene().removeItem(self.cut_line_rubber_band)

        del self._toolbar

    def on_action_divide_triggered(self):
        layers = self._iface.layerTreeView().selectedLayers()

        if len(layers) == 0:
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'No selected layer!')
            return

        active_layer = self._iface.activeLayer()

        if active_layer is None:
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'No active layer!')
            return

        if active_layer.__class__.__name__ != 'QgsVectorLayer':
            QMessageBox.information(self._iface.mainWindow(), 'Info', 'Selected layer is not Vector layer!')
            return

        geom_type = active_layer.geometryType()

        if geom_type is not None:
            if geom_type != QgsWkbTypes.PolygonGeometry:
                QMessageBox.information(self._iface.mainWindow(), 'Info', 'Selected layer is not Polygon layer!')
                return

        self.map_canvas.setMapTool(self.reference_line_selector)

    def on_action_reference_line_selector_triggered(self):
        self.map_canvas.setMapTool(self.reference_line_selector)

    def hide_all_marker_rubber_band(self):
        self.reference_line_marker.hide()
        self.line_rubber_band.hide()

    def on_reference_line_selector_canvas_clicked(self, point, button):
        """
        :type: QgsPointXY
        :param: point

        """

        active_layer = self._iface.activeLayer()

        point_geometry = QgsGeometry.fromPointXY(point)
        nearest_polygon_feature = find_nearest_feature(point_geometry, active_layer)
        dist = nearest_polygon_feature.geometry().distance(point_geometry)
        tolerance = 1

        log_message("dist " + str(dist))

        if dist > tolerance:
            log_message("failed to find nearest_polygon_feature")
            self.reference_line_marker.hide()
            return

        # get start and end of nearest edge of the nearest polygon

        nearest_edge = find_nearest_edge_from_multi_polygon(point_geometry, nearest_polygon_feature.geometry())

        start_point = nearest_edge[0]
        end_point = nearest_edge[1]

        # get mid point of the nearest edge
        mid_x = (start_point.x() + end_point.x()) / 2
        mid_y = (start_point.y() + end_point.y()) / 2

        self.reference_line_marker.setCenter(QgsPointXY(mid_x, mid_y))

        self.reference_line_marker.show()
        nearest_edge = QgsGeometry.fromPolylineXY([start_point, end_point])

        self.line_rubber_band.setToGeometry(nearest_edge)
        self.line_rubber_band.show()

        distance, ok = QInputDialog.getDouble(self._iface.mainWindow(), "Input", "enter a distance")

        if not ok:
            self.hide_all_marker_rubber_band()
            return

        # get translated edge
        offset_edge = nearest_edge.offsetCurve(distance, 8, QgsGeometry.JoinStyleRound, 5.0)

        offset_edge_start = offset_edge.asPolyline()[0]
        offset_edge_end = offset_edge.asPolyline()[1]

        bounding_box = nearest_polygon_feature.geometry().boundingBox()

        extend_length = max(bounding_box.width(), bounding_box.height())

        edge_length = distance_between_two_points(start_point, end_point)

        # get mid point of the offset_edge
        mid_x = (offset_edge_start.x() + offset_edge_end.x()) / 2
        mid_y = (offset_edge_start.y() + offset_edge_end.y()) / 2

        # extend offset_edge by edge_length / 2
        new_start_x = offset_edge_start.x() + (offset_edge_start.x() - mid_x) / (edge_length / 2) * extend_length
        new_start_y = offset_edge_start.y() + (offset_edge_start.y() - mid_y) / (edge_length / 2) * extend_length

        new_end_x = offset_edge_end.x() + (offset_edge_end.x() - mid_x) / (edge_length / 2) * extend_length
        new_end_y = offset_edge_end.y() + (offset_edge_end.y() - mid_y) / (edge_length / 2) * extend_length

        cut_line = QgsGeometry.fromPolylineXY([QgsPointXY(new_start_x, new_start_y), QgsPointXY(new_end_x, new_end_y)])

        if self._show_cut_line:
            self.cut_line_rubber_band.setToGeometry(cut_line)
            self.cut_line_rubber_band.show()

        # split polygon by line

        geom = nearest_polygon_feature.geometry()
        success, split_geometry_list, topo = geom.splitGeometry(cut_line.asPolyline(), False)

        if success != QgsGeometry.OperationResult.Success:
            self.hide_all_marker_rubber_band()
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "failed to split")
            return

        if len(split_geometry_list) != 1:
            self.hide_all_marker_rubber_band()
            QMessageBox.critical(self._iface.mainWindow(), 'Error', "failed to split")
            return

        active_layer.startEditing()

        active_layer.changeGeometry(nearest_polygon_feature.id(), geom)

        feature = QgsVectorLayerUtils.createFeature(active_layer)

        feature.setGeometry(split_geometry_list[0])
        active_layer.addFeature(feature)

        active_layer.commitChanges()


if __name__ == '__main__':
    pass
