package com.Sts.Actions.Crossplot;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
  * <p>Description: </p>
  * <p>Copyright: Copyright (c) 2002</p>
  * <p>Company: S2S Systems LLC</p>
  *
  * @author TJLasseter
  * @version 1.1
  */

 public class StsXPolygonAction extends StsAction
 {
     StsCrossplotClass crossplotClass = null;
     StsCrossplot crossplot;
     StsSeismicVolume[] seismicVolumes;
     byte[][] planesData;
     StsSeismicCursorSection seismicCursorSection;
     StsColorItemComboBoxFieldBean colorItemListBean = null;
     int nRows, nCols;
     StsColor currentColor = null;
     String currentName = null;

     int nPicks = 0;
     StsPoint firstPick = null;
     StsPoint pick;
     StsToggleButton polygonOrPointToggleButton = null;

     static public StsXPolygon selectedPolygon = null;
     static public StsXPolygon pressedPolygon = null;
     static public StsXPolygon releasedPolygon = null;

     static public final int NONE = 0;

     static public final int MODE_CREATE = 1;
     static public final int MODE_EDIT = 2;
     static public final String[] modeLabels = new String[]{"NONE", "MODE_CREATE", "MODE_EDIT"};

     static public final int TYPE_EDGE = 1;
     static public final int TYPE_VERTEX = 2;
     static public final String[] typeLabels = new String[]{"NONE", "TYPE_EDGE", "TYPE_VERTEX"};

     static public int mode = NONE;
     static public int type = NONE;

     static public boolean onFirstPick = false;
     static public int pickIndex = -1;
     static public double[] savedPoint;

     static StsXPolygon[] seismicPolygons = null;
     static int dirNo = -1;
     static float dirCoordinate = StsParameters.nullValue;
     static StsCrossplotPoint[] seismicCrossplotPoints = null;
     static StsCrossplotPoint[] seismicCrossplotPolygonPoints = null;
     static boolean polygonAxesFlipped = false;

     static boolean debug = false;

     public StsXPolygonAction(StsActionManager actionManager)
     {
         super(actionManager);
     }

     public boolean start()
     {
         try
         {
             crossplotClass = StsCrossplot.getCrossplotClass();
             crossplot = crossplotClass.getCurrentCrossplot();
             if (crossplot == null)
             {
                 StsException.systemError("StsXPolygonAction.constructor() failed: crossplot is null.");
                 deselectButton();
                 return false;
             }

             seismicVolumes = crossplot.getVolumes();
             if (seismicVolumes == null || seismicVolumes.length < 2)
             {
                 StsException.systemError("StsXPolygonAction.constructor() failed: number of seismic volumes < 2.");
                 deselectButton();
                 return false;
             }
             initializeViewCursorPicks();
             initializePicks();

             colorItemListBean = (StsColorItemComboBoxFieldBean) model.win3d.getToolbarComponentNamed(StsCrossplotToolbar.NAME, StsCrossplotToolbar.EDIT_COLORS);
             colorItemListBean.setEditable(true);
             colorItemListBean.initialize(crossplot, "typeItem", null, crossplot.getTypeLibrary().getColorListItems());
             initializeCurrentColor();

             polygonOrPointToggleButton = (StsToggleButton) model.win3d.getToolbarComponentNamed(StsCrossplotToolbar.NAME, StsCrossplotToolbar.SELECT_POINT_OR_POLYGON);
             //            viewChanged(glPanel3d.getCurrentView());

             logMessage("Pick polygon points on crossplot or crossplot points on seismic.");
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsXPolygonAction.start() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     private void initializeCurrentColor()
     {
         StsColorListItem item = colorItemListBean.getItemAt(0);
         StsType type = (StsType) item.getObject();
         currentColor = type.getStsColor();
     }

     private void initializeViewCursorPicks()
     {
         seismicVolumes = crossplot.getVolumes();
         StsCursor3d cursor3d = model.win3d.getCursor3d();
         dirNo = cursor3d.getCurrentDirNo();
         dirCoordinate = cursor3d.getCurrentDirCoordinate();
         //        nPlane = glPanel3d.getCursor3d().getCurrentPlaneIndex();
         nRows = seismicVolumes[0].getNCursorRows(dirNo);
         nCols = seismicVolumes[0].getNCursorCols(dirNo);

         planesData = new byte[2][];
         for (int n = 0; n < 2; n++)
             planesData[n] = seismicVolumes[n].readBytePlaneData(dirNo, dirCoordinate);
     }

     private void initializePicks()
     {
         //        crossplot.setCurrentCrossplotColors();
         initializePolygon();
     }

     private void initializePolygon()
     {
         selectedPolygon = null;
         mode = NONE;
         type = NONE;
         onFirstPick = false;
         pickIndex = -1;
     }

     /** user has toggled button down to start action. It has failed, so we must toggle it up. */
     private void deselectButton()
     {
         model.win3d.selectToolbarItem(StsCrossplotToolbar.NAME, StsCrossplotToolbar.CREATE_POLYGON, false);
     }

     public void viewChanged(Object newView)
     {
         polygonOrPointToggleButton.setEnabled(newView instanceof StsViewCursor);
     }

     public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
     {
         StsView currentView = ((StsGLPanel3d) glPanel).getView();
         if (currentView instanceof StsViewXP)
             return performMousePolygonAction(mouse, (StsViewXP) currentView, glPanel);
         else if (currentView instanceof StsViewCursor)
         {
             if (polygonOrPointToggleButton.isSelected()) // if down (selected) do points, otherwise do polygons
             {
                 if (selectedPolygon != null)
                 {
                     if (!selectedPolygon.isClosed())
                     {
                         deleteSeismicPolygon(selectedPolygon);
                         initializePolygon();
                     }
                 }
                 return performMouseAction(mouse, (StsViewCursor) currentView, (StsGLPanel3d) glPanel);
             }
             else
                 return performMousePolygonAction(mouse, (StsViewCursor) currentView, (StsGLPanel3d) glPanel);
         }
         else if (currentView instanceof StsView3d)
         {
             //            if(polygonOrPointToggleButton.isSelected())
             return performMouseAction(mouse, (StsView3d) currentView, (StsGLPanel3d) glPanel);
             //            else
             //                return performMousePolygonAction(mouse, (StsView3d)currentView);
         }
         else // a pick in the 3d window does nothing
             return true;
     }

     /** mouse action for 3d crossplot window: return false ONLY if action fails completely. */
     private boolean performMousePolygonAction(StsMouse mouse, StsViewXP viewXP, StsGLPanel glPanel)
     {
         if (crossplot != crossplotClass.getCurrentCrossplot())
         {
             if (selectedPolygon != null) selectedPolygon.close();
             crossplot.clearCrossplotTextureDisplays();
             crossplot = crossplotClass.getCurrentCrossplot();
             model.viewObjectRepaint(this, crossplot);
             return true;
         }

         //        viewXP.insetViewPort();
         pick = viewXP.computePickPoint(mouse);
         boolean status = handlePolygonAction(mouse, viewXP, glPanel);
         //        viewXP.resetViewPort();

         return status;
     }

     private boolean handlePolygonAction(StsMouse mouse, StsViewXP viewXP, StsGLPanel glPanel)
     {
         int buttonState;

         //        view2d = (StsView2d) glPanel3d.getCurrentView();
         //        view2d.insetViewPort();

         if (mouse.getCurrentButton() == StsMouse.LEFT)
         {
             buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
             onFirstPick = false;

             if (debug)
             {
                 mouse.printState();
                 printMode();
                 printType();
             }

             if (buttonState == StsMouse.PRESSED)
             {
                 pressedPolygon = (StsXPolygon) StsJOGLPick.pickClass3d(glPanel, crossplot.getPolygons(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_ALL);

                 if (pressedPolygon != null && pressedPolygon != selectedPolygon)
                 {
                     getPolygonLibraryType(pressedPolygon);
                 }

                 if (mode == NONE)
                 {
                     if (pressedPolygon == null)
                     {
                         getCrossplotLibraryType();
                         if (currentColor == null)
                         {
                             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
                             return true;
                         }
                         setMode(MODE_CREATE);
                         firstPick = pick;
                         selectedPolygon = null;
                         addCrossplotPolygonPoint();
                     }
                     else
                     {
                         setMode(MODE_EDIT);
                         selectedPolygon = pressedPolygon;
                         setType(TYPE_EDGE);
                     }
                 }
                 else if (mode == MODE_CREATE)
                 {
                     addCrossplotPolygonPoint();
                     onFirstPick = nPicks >= 3 && viewXP.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (pressedPolygon == null)
                     {
                         selectedPolygon = null;
                         pickIndex = -1;
                         setType(NONE);
                         setMode(NONE);
                     }
                     else if (pressedPolygon != selectedPolygon)  // selectedPolygon deselected, pressed one selected
                     {
                         selectedPolygon = pressedPolygon;
                         pickIndex = -1;
                         setType(NONE);
                     }
                     else // pick on previously selected polygon
                     {
                         //                        int previousType = type;
                         //                        int previousPickIndex = pickIndex;

                         StsPickItem[] pickItems = StsJOGLPick.pickItems;
                         type = pickItems[0].names[1];
                         pickIndex = pickItems[0].names[2];
                         if (type == TYPE_EDGE) // otherwise if pick on edge, we are inserting a vertex
                         {
                             if (debug) System.out.println("Edge " + pickIndex + " picked.");
                             pickIndex = selectedPolygon.insertPoint(pickIndex, pick, debug);
                             setType(TYPE_VERTEX);
                         }
                         else
                         {
                             if (debug) System.out.println("Vertex " + pickIndex + " picked.");
                             double[] pickedPoint = selectedPolygon.getPoint(pickIndex);
                             savedPoint = new double[pickedPoint.length];
                             System.arraycopy(pickedPoint, 0, savedPoint, 0, pickedPoint.length);
                         }
                         pressedPolygon = null;
                     }

                 }
             }
             // Rubber band from last or between two on either side
             else if (buttonState == StsMouse.DRAGGED)
             {
                 //                pressedPolygon = null;
                 logMessage("X: " + pick.v[0] + " Y: " + pick.v[1]);
                 if (mode == MODE_EDIT && type == TYPE_VERTEX)
                 {
                     selectedPolygon.movePoint(pickIndex, pick, debug);
                     crossplot.polygonsChanged();
                     model.viewObjectChanged(this, crossplot);
                 }
                 else if (mode == MODE_CREATE)
                 {
                     selectedPolygon.moveLastPoint(pick, debug);
                     onFirstPick = nPicks >= 3 && viewXP.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (debug) System.out.println("Selected point being dragged.");
                 }
                 //                viewXP.glPanel3d.repaint();
             }
             else if (buttonState == StsMouse.RELEASED) // add a point unless near firstPoint && nPicks > 3
             {
                 if (mode == MODE_CREATE)
                 {
                     onFirstPick = nPicks >= 3 && viewXP.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (onFirstPick)
                     {
                         selectedPolygon.deleteLastPoint();
                         endPolygon();
                         setMode(NONE);
                         if (debug) System.out.println("Last point on first point: polygon joined.");
                     }
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (type == TYPE_VERTEX)
                     {
                         if (!viewXP.isInsideInsetViewPort(mouse.mousePoint.x, mouse.mousePoint.y))
                             selectedPolygon.movePoint(pickIndex, savedPoint, false);
                         model.instanceChange(selectedPolygon, "point moved");
                         crossplot.clearCrossplotTextureDisplays();
                         crossplot.processPolygons();
                     }
                 }
                 model.viewObjectChanged(this, crossplot);
             }
         }
         else if (mouse.getCurrentButton() == StsMouse.VIEW)
         {
             buttonState = mouse.getButtonStateCheckClear(StsMouse.VIEW);

             // Delete Vertex
             if (buttonState == StsMouse.RELEASED)
             {
                 // Delete the entire polygon
                 if (nPicks < 3)           // Delete Entire Polygon if only 3 points left
                 {

                 }
                 else                      // Delete vertex from polygon
                 {

                 }
             }
         }
         viewXP.resetViewPort();
         model.viewObjectRepaint(this, crossplot);
         return true;
     }

     private StsType getPolygonLibraryType(StsXPolygon pressedPolygon)
     {
         StsType libraryType = pressedPolygon.getStsType();
         if (libraryType == null) return null;
         colorItemListBean.setSelectedIndex(crossplot.getTypeLibrary().getTypes().getIndex(libraryType));
         currentColor = libraryType.getStsColor();
         currentName = libraryType.getName();
         return libraryType;
     }

     private StsType getCrossplotLibraryType()
     {
         StsTypeLibrary typeLibrary = crossplot.getTypeLibrary();
         if (typeLibrary == null) return null;
         StsType libraryType = typeLibrary.getCurrentType();
         if (libraryType == null) return null;
         currentColor = libraryType.getStsColor();
         currentName = libraryType.getName();
         return typeLibrary.getCurrentType();
     }

     private void setMode(int mode)
     {
         this.mode = mode;
         if (debug) System.out.println("Mode: " + modeLabels[mode]);
     }

     private void setType(int type)
     {
         this.type = type;
         if (debug) System.out.println("Type: " + typeLabels[type]);
     }

     private void printMode() { System.out.println("Mode: " + modeLabels[mode]); }

     private void printType() { System.out.println("Type: " + typeLabels[type]); }

     public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLPanel glPanel)
     {
         if (mouse.isButtonDown()) return false;
         if (selectedPolygon == null) return false;

         StsGLPanel3d glPanel3d = (StsGLPanel3d) glPanel;
         char keyReleased = e.getKeyChar();
         if (keyReleased == 'D' || keyReleased == 'd') // vertex or polygon is being deleted
         {
             StsView currentView = glPanel3d.getView();

             int nPoints = selectedPolygon.getNPoints();
             if (type == TYPE_VERTEX && nPoints > 3)
             {
                 if (currentView instanceof StsViewXP)
                 {
                     crossplot.deletePolygonPoint(selectedPolygon, pickIndex);
                     model.instanceChange(selectedPolygon, "crossplot polygon point deleted");
                     crossplot.clearCrossplotTextureDisplays();
                     crossplot.processPolygons();
                     model.viewObjectChangedAndRepaint(this, crossplot);
                     return true;
                 }
                 else
                 {
                     deletePolygonPoint(selectedPolygon, pickIndex);
                     processSeismicPolygons();
                     model.viewObjectChanged(this, crossplot);
                 }
                 if (pickIndex == selectedPolygon.getNPoints()) pickIndex = 0;
             }
             else if (type == TYPE_EDGE || nPoints <= 3)
             {
                 if (currentView instanceof StsViewXP)
                 {
                     crossplot.deletePolygon(selectedPolygon);
                     endAndStartTransaction();
                     crossplot.clearCrossplotTextureDisplays();
                     crossplot.processPolygons();
                     model.viewObjectChanged(this, crossplot);
                 }
                 else
                 {
                     deleteSeismicPolygon(selectedPolygon);
                     model.viewObjectChanged(this, crossplot);
                 }

                 pickIndex = -1;
                 type = NONE;
                 mode = NONE;
             }
             else
             {
                 mode = NONE;
                 return true;
             }
         }
         return false;
     }

     private void repaintCrossplotViews()
     {
         model.repaintViews(StsView3d.class);
         model.repaintViews(StsViewCursor.class);
         model.repaintViews(StsViewXP.class);
     }

     private void crossplotChangedRepaint()
     {
         model.viewObjectChangedAndRepaint(this, crossplot);
     }

     private void endAndStartTransaction()
     {
         model.commit();
         model.getCreateCurrentTransaction("StsXPolygonAction");
     }

     /** mouse action for 2d seismic window: return false ONLY if action fails completely. */
     private boolean performMousePolygonAction(StsMouse mouse, StsViewCursor viewCursor, StsGLPanel glPanel)
     {
         checkCurrentSeismicPolygons();
         return handlePolygonAction(mouse, viewCursor, glPanel);
     }

     private void checkCurrentSeismicPolygons()
     {
         StsCursor3d cursor3d = model.win3d.getCursor3d();
         if (dirNo != cursor3d.getCurrentDirNo() || dirCoordinate != cursor3d.getCurrentDirCoordinate())
         {
             this.clearSeismicPoints();
             initializePolygon();
             StsSeismicCursorSection seismicCursorSection = (StsSeismicCursorSection) cursor3d.getDisplayableSection(dirNo, StsSeismicCursorSection.class);
             dirNo = seismicCursorSection.dirNo;
             dirCoordinate = seismicCursorSection.dirCoordinate;
             //        nPlane = seismicCursorSection.nPlane;
         }
     }

     private boolean handlePolygonAction(StsMouse mouse, StsViewCursor viewCursor, StsGLPanel glPanel)
     {
         int buttonState;

         pick = viewCursor.computePickPoint(mouse);
         if (pick == null) return true;

         polygonAxesFlipped = viewCursor.axesFlipped;

         if (mouse.getCurrentButton() == StsMouse.LEFT)
         {
             buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
             //            if(mode == NONE) polygons = seismicCursorSection.getSeismicPolygons();
             onFirstPick = false;

             if (debug)
             {
                 mouse.printState();
                 printMode();
                 printType();
                 pick.print();
             }

             if (buttonState == StsMouse.PRESSED)
             {
                 pressedPolygon = (StsXPolygon) StsJOGLPick.pickClass3d(glPanel, seismicPolygons, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_ALL);

                 if (pressedPolygon != null && pressedPolygon != selectedPolygon)
                 {
                     getPolygonLibraryType(pressedPolygon);
                 }

                 if (mode == NONE)
                 {
                     if (pressedPolygon == null)
                     {
                         getCrossplotLibraryType();
                         if (currentColor == null)
                         {
                             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
                             return true;
                         }
                         setMode(MODE_CREATE);
                         firstPick = pick;
                         selectedPolygon = null;
                         addSeismicPolygonPoint();
                     }
                     else
                     {
                         setMode(MODE_EDIT);
                         selectedPolygon = pressedPolygon;
                         setType(TYPE_EDGE);
                     }
                 }
                 else if (mode == MODE_CREATE)
                 {
                     addSeismicPolygonPoint();
                     onFirstPick = nPicks >= 3 && viewCursor.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (pressedPolygon == null)
                     {
                         selectedPolygon = null;
                         pickIndex = -1;
                         setType(NONE);
                         setMode(NONE);
                     }
                     else if (pressedPolygon != selectedPolygon)  // selectedPolygon deselected, pressed one selected
                     {
                         selectedPolygon = pressedPolygon;
                         pickIndex = -1;
                         setType(NONE);
                     }
                     else // pick on previously selected polygon
                     {
                         //                        int previousType = type;
                         //                        int previousPickIndex = pickIndex;

                         StsPickItem[] pickItems = StsJOGLPick.pickItems;
                         type = pickItems[0].names[1];
                         pickIndex = pickItems[0].names[2];
                         if (type == TYPE_EDGE) // otherwise if pick on edge, we are inserting a vertex
                         {
                             if (debug) System.out.println("Edge " + pickIndex + " picked.");
                             pickIndex = selectedPolygon.insertPoint(pickIndex, pick, debug);
                             setType(TYPE_VERTEX);
                         }
                         else if (debug) System.out.println("Vertex " + pickIndex + " picked.");

                         pressedPolygon = null;
                     }

                 }
             }
             // Rubber band from last or between two on either side
             else if (buttonState == StsMouse.DRAGGED)
             {
                 //                pressedPolygon = null;
                 logMessage("X: " + pick.v[0] + " Y: " + pick.v[1]);
                 if (mode == MODE_EDIT && type == TYPE_VERTEX)
                 {
                     selectedPolygon.movePoint(pickIndex, pick, debug);
                     processSeismicPolygons();
                 }
                 else if (mode == MODE_CREATE)
                 {
                     selectedPolygon.moveLastPoint(pick, debug);
                     onFirstPick = nPicks >= 3 && viewCursor.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (debug) System.out.println("Selected point being dragged.");
                 }
                 viewCursor.glPanel3d.repaint();
                 model.viewObjectRepaint(this, crossplot);
             }
             else if (buttonState == StsMouse.RELEASED) // add a point unless near firstPoint && nPicks > 3
             {
                 if (mode == MODE_CREATE)
                 {
                     onFirstPick = nPicks >= 3 && viewCursor.mousePicksNear(pick, firstPick, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (onFirstPick)
                     {
                         selectedPolygon.deleteLastPoint();
                         endSeismicPolygon();
                         setMode(NONE);
                         if (debug) System.out.println("Last point on first point: polygon joined.");
                     }
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (type == TYPE_VERTEX)
                     {
                         //                        crossplot.clearCrossplotDataDisplays();
                         //                        crossplot.processPolygons();
                     }
                 }
                 crossplotChangedRepaint();
             }
         }
         else if (mouse.getCurrentButton() == StsMouse.VIEW)
         {
             buttonState = mouse.getButtonStateCheckClear(StsMouse.VIEW);

             // Delete Vertex
             if (buttonState == StsMouse.RELEASED)
             {
                 // Delete the entire polygon
                 if (nPicks < 3)           // Delete Entire Polygon if only 3 points left
                 {

                 }
                 else                      // Delete vertex from polygon
                 {

                 }
             }
         }
         viewCursor.resetViewPort();
         return true;
     }

     private boolean performMouseAction(StsMouse mouse, StsViewCursor viewCursor, StsGLPanel3d glPanel3d)
     {
         StsCursorPoint cursorPoint = null;

         if (dirNo == StsParameters.NO_MATCH) return false;

         StsCursor3d cursor3d = model.win3d.getCursor3d();
         StsSeismicCursorSection seismicCursorSection = (StsSeismicCursorSection)cursor3d.getDisplayableSection(dirNo, StsCrossplot.class);
         if (seismicCursorSection == null) return false;

         if (this.seismicCursorSection != seismicCursorSection)
         {
             initializeViewCursorPicks();
         }

         if (mouse.getCurrentButton() != StsMouse.LEFT) return false;

         getCrossplotLibraryType();
         if (currentColor == null)
         {
             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
             return true;
         }

         int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

         pick = viewCursor.computePickPoint3d(mouse);
         if (pick == null) return true;
         float[] dataRowColF = seismicCursorSection.getDataRowColF(pick);
         if (dataRowColF == null) return true;
         int nDir = glPanel3d.getCursor3d().getCurrentDirNo();
         cursorPoint = new StsCursorPoint(nDir, dataRowColF, pick);
         if (cursorPoint == null) return true;

         if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
         {
             String readout = seismicCursorSection.logReadout2d(dataRowColF);
             logMessage(readout);
         }
         else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
         {
             if (cursorPoint != null)
             {
                 int row = Math.round(dataRowColF[0]);
                 int col = Math.round(dataRowColF[1]);
                 int index = row * nCols + col;
                 byte colValue = planesData[0][index];
                 byte rowValue = planesData[1][index];
                 addCrossplotPoint(cursorPoint, row, col, currentColor, colValue, rowValue);

                 String readout = seismicCursorSection.logReadout2d(dataRowColF);
                 logMessage("Added crossplot point: " + " " + readout);
                 //                model.win3d.win3dDisplay();
                 crossplotChangedRepaint();
             }
         }
         return true;
     }
/*
     private StsSeismicCursorSection getSeismicCursorSection(StsGLPanel3d glPanel3d)
     {
         StsCursor3dTexture[] cursor3dTextures = glPanel3d.getDisplayableSections(dirNo);
         for (StsCursor3dTexture cursor3dTexture : cursor3dTextures)
         {
             if (cursor3dTexture instanceof StsSeismicCursorSection)
                 return (StsSeismicCursorSection) cursor3dTexture;
         }
         return null;
     }
*/
     //TODO DEBUG:  confusion in StsCursorPoint and here on data row/col, cursor row/col, volume row/col
     private void addCrossplotPoint(StsCursorPoint cursorPoint, int dataRow, int dataCol, StsColor currentColor, byte colValue, byte rowValue)
     {
         StsCrossplotPoint crossplotPoint = new StsCrossplotPoint(cursorPoint, dataRow, dataCol, currentColor, colValue, rowValue);
         seismicCrossplotPoints = (StsCrossplotPoint[]) StsMath.arrayAddElement(seismicCrossplotPoints, crossplotPoint);
     }

     /** mouse action for 3d seismic window: return false ONLY if action fails completely. */
/*
     private boolean performMousePolygonAction(StsMouse mouse, StsView3d view3d, StsGLPanel glPanel)
     {
         checkCurrentSeismicPolygons();
         return handlePolygonAction(mouse, view3d, (StsGLPanel3d) glPanel);
     }

     private boolean handlePolygonAction(StsMouse mouse, StsView3d view3d, StsGLPanel3d glPanel3d)
     {
         int buttonState;

         if (mouse.getCurrentButton() != StsMouse.LEFT) return false;

         getCrossplotLibraryType();

         if (currentColor == null)
         {
             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
             return true;
         }

         int leftButtonState = mouse.getButtonState(StsMouse.LEFT);

         int nDir = glPanel3d.getCursor3d().getCurrentDirNo();
         pick = glPanel3d.getCursor3d().getPointInCursorPlane(glPanel3d, nDir, mouse);
         if (pick == null) return true;

         if (mouse.getCurrentButton() == StsMouse.LEFT)
         {
             buttonState = mouse.getButtonState(StsMouse.LEFT);
             //            if(mode == NONE) polygons = seismicCursorSection.getSeismicPolygons();
             onFirstPick = false;

             if (debug)
             {
                 mouse.printState();
                 printMode();
                 printType();
             }

             if (buttonState == StsMouse.PRESSED)
             {
                 pressedPolygon = (StsXPolygon) StsJOGLPick.pickClass3d(glPanel3d, seismicPolygons, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_ALL);

                 if (pressedPolygon != null && pressedPolygon != selectedPolygon)
                 {
                     getPolygonLibraryType(pressedPolygon);
                 }

                 if (mode == NONE)
                 {
                     if (pressedPolygon == null)
                     {
                         getCrossplotLibraryType();
                         if (currentColor == null)
                         {
                             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
                             return true;
                         }
                         setMode(MODE_CREATE);
                         firstPick = pick;
                         selectedPolygon = null;
                         addSeismicPolygonPoint();
                     }
                     else
                     {
                         setMode(MODE_EDIT);
                         selectedPolygon = pressedPolygon;
                         setType(TYPE_EDGE);
                     }
                 }
                 else if (mode == MODE_CREATE)
                 {
                     addSeismicPolygonPoint();
                     onFirstPick = nPicks >= 3 && StsJOGLPick.mousePicksPoint(mouse.getMousePoint(), firstPick, glPanel3d, StsJOGLPick.PICKSIZE_MEDIUM);
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (pressedPolygon == null)
                     {
                         selectedPolygon = null;
                         pickIndex = -1;
                         setType(NONE);
                         setMode(NONE);
                     }
                     else if (pressedPolygon != selectedPolygon)  // selectedPolygon deselected, pressed one selected
                     {
                         selectedPolygon = pressedPolygon;
                         pickIndex = -1;
                         setType(NONE);
                     }
                     else // pick on previously selected polygon
                     {
                         //                        int previousType = type;
                         //                        int previousPickIndex = pickIndex;

                         StsPickItem[] pickItems = StsJOGLPick.pickItems;
                         type = pickItems[0].names[1];
                         pickIndex = pickItems[0].names[2];
                         if (type == TYPE_EDGE) // otherwise if pick on edge, we are inserting a vertex
                         {
                             if (debug) System.out.println("Edge " + pickIndex + " picked.");
                             pickIndex = selectedPolygon.insertPoint(pickIndex, pick, debug);
                             setType(TYPE_VERTEX);
                         }
                         else if (debug) System.out.println("Vertex " + pickIndex + " picked.");

                         pressedPolygon = null;
                     }

                 }
             }
             // Rubber band from last or between two on either side
             else if (buttonState == StsMouse.DRAGGED)
             {
                 //                pressedPolygon = null;
                 logMessage("X: " + pick.v[0] + " Y: " + pick.v[1]);
                 if (mode == MODE_EDIT && type == TYPE_VERTEX)
                 {
                     selectedPolygon.movePoint(pickIndex, pick, debug);
                     processSeismicPolygons();
                 }
                 else if (mode == MODE_CREATE)
                 {
                     selectedPolygon.moveLastPoint(pick, debug);
                     onFirstPick = nPicks >= 3 && StsJOGLPick.mousePicksPoint(mouse.getMousePoint(), firstPick, glPanel3d, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (debug) System.out.println("Selected point being dragged.");
                 }
                 view3d.glPanel3d.repaint();
             }
             else if (buttonState == StsMouse.RELEASED) // add a point unless near firstPoint && nPicks > 3
             {
                 if (mode == MODE_CREATE)
                 {
                     onFirstPick = nPicks >= 3 && StsJOGLPick.mousePicksPoint(mouse.getMousePoint(), firstPick, glPanel3d, StsJOGLPick.PICKSIZE_MEDIUM);
                     if (onFirstPick)
                     {
                         selectedPolygon.deleteLastPoint();
                         endSeismicPolygon();
                         setMode(NONE);
                         if (debug) System.out.println("Last point on first point: polygon joined.");
                     }
                 }
                 else if (mode == MODE_EDIT)
                 {
                     if (type == TYPE_VERTEX)
                     {
                         //                        crossplot.clearCrossplotDataDisplays();
                         //                        crossplot.processPolygons();
                     }
                 }
                 mouse.clearButtonState(StsMouse.LEFT);
                 crossplotChangedRepaint();
             }
         }
         else if (mouse.getCurrentButton() == StsMouse.VIEW)
         {
             buttonState = mouse.getButtonState(StsMouse.VIEW);

             // Delete Vertex
             if (buttonState == StsMouse.RELEASED)
             {
                 // Delete the entire polygon
                 if (nPicks < 3)           // Delete Entire Polygon if only 3 points left
                 {

                 }
                 else                      // Delete vertex from polygon
                 {

                 }
                 mouse.clearButtonState(StsMouse.LEFT);
             }
         }
         return true;
     }
*/
     private boolean performMouseAction(StsMouse mouse, StsView3d view3d, StsGLPanel3d glPanel3d)
     {
         StsCursorPoint cursorPoint = null;

         if (dirNo == StsParameters.NO_MATCH) return false;
         StsCursor3d cursor3d = model.win3d.getCursor3d();
         StsSeismicCursorSection seismicCursorSection = (StsSeismicCursorSection)cursor3d.getDisplayableSection(dirNo, StsCrossplot.class);
         if (seismicCursorSection == null) return false;

         if (this.seismicCursorSection != seismicCursorSection)
         {
             this.seismicCursorSection = seismicCursorSection;
             //            seismicCursorSection.clearCrossplotPoints();
             clearSeismicPoints();
             initializeViewCursorPicks();
         }

         if (mouse.getCurrentButton() != StsMouse.LEFT) return false;

         getCrossplotLibraryType();

         if (currentColor == null)
         {
             new StsMessage(model.win3d, StsMessage.INFO, "No polygon color selected; select and try again.");
             return true;
         }

         int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

         int nDir = cursor3d.getCurrentDirNo();
         cursorPoint = cursor3d.getNearestPointInCursorPlane(glPanel3d, mouse, seismicVolumes[0]);
         if (cursorPoint == null) return true;
         if (cursorPoint.dirNo != nDir)
         {
             new StsMessage(model.win3d, StsMessage.WARNING,
                 "Picked point is not on current section which is the " + StsCursor3d.coorLabels[nDir] + " cursor.");
             return true;
         }

         pick = cursor3d.getPointInCursorPlane(glPanel3d, nDir, mouse);
         if (pick == null) return true;
         float[] dataRowColF = seismicCursorSection.getDataRowColF(pick);
         if (dataRowColF == null) return true;
         int dataRow = Math.round(dataRowColF[0]);
         int dataCol = Math.round(dataRowColF[1]);
         cursorPoint = new StsCursorPoint(nDir, pick);

         if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
             cursor3d.logReadout(glPanel3d, cursorPoint);
         else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
         {
             if (cursorPoint != null)
             {
                 int row = Math.round(cursorPoint.rowNum);
                 int col = Math.round(cursorPoint.colNum);
                 int index = row * nCols + col;
                 byte colValue = planesData[0][index];
                 byte rowValue = planesData[1][index];
                 addCrossplotPoint(cursorPoint, dataRow, dataCol, currentColor, colValue, rowValue);
                 cursor3d.logReadout(glPanel3d, "Added crossplot point: ", cursorPoint);
                 //                model.win3d.win3dDisplay();
                 crossplotChangedRepaint();
             }
         }
         return true;
     }

     public boolean end()
     {
         if (seismicCursorSection != null)
         {
             //            seismicCursorSection.clearCrossplotPoints();
             clearSeismicPoints();
         }
         crossplot.checkPolygons();
         statusArea.textOnly();
         return true;
     }

     private void clearSelection()
     {
         selectedPolygon = null;
         //        mode = NONE;
         pickIndex = -1;
         type = NONE;
     }

     private void addCrossplotPolygonPoint()
     {
         if (selectedPolygon == null)
         {
             StsType polygonType = crossplot.getPolygonType(currentName, currentColor);
             selectedPolygon = new StsXPolygon(polygonType);
             if (debug) System.out.println("new polygon created.");
             crossplot.addPolygon(selectedPolygon);
         }
         selectedPolygon.addPoint(pick, debug);
         nPicks++;
         logMessage("Added point to current polygon.");
     }

     private void addSeismicPolygonPoint()
     {
         if (selectedPolygon == null)
         {
             StsType polygonType = crossplot.getPolygonType(currentName, currentColor);
             selectedPolygon = new StsXPolygon(polygonType, false); // polygons are non-persistent
             if (debug) System.out.println("new polygon created.");
             addSeismicPolygon(selectedPolygon);
         }
         selectedPolygon.addPoint(pick, debug);
         nPicks++;
         logMessage("Added point to current polygon.");
     }

     private void endPolygon()
     {
         selectedPolygon.close();
         endAndStartTransaction();
         crossplot.clearCrossplotTextureDisplays();
         crossplot.processPolygons();
         selectedPolygon = null;
         nPicks = 0;
     }

     private void endSeismicPolygon()
     {
         selectedPolygon.close();
         selectedPolygon = null;
         nPicks = 0;
         processSeismicPolygons();
     }

     public boolean processSeismicPolygons()
     {
         StsGLOffscreenPolygon offscreen;

         try
         {
             seismicCrossplotPolygonPoints = null;

             if (seismicPolygons == null) return false;
             int nPolygons = seismicPolygons.length;
             if (nPolygons == 0) return false;
             if (crossplot == null) return false;

             int[] nRowCols = seismicVolumes[0].getCursorDataNRowCols(dirNo);
             int nRows = nRowCols[0];
             int nCols = nRowCols[1];

             float[][] cursorRange = seismicVolumes[0].getCursorAxisRange(dirNo);
             offscreen = new StsGLOffscreenPolygon(nCols, nRows, cursorRange[0][0], cursorRange[0][1], cursorRange[1][0], cursorRange[1][1], false);

             int[] coorIndexes = StsSeismicVolume.getCursor2dCoorDataIndexes(dirNo, polygonAxesFlipped);

             offscreen.setPolygons(seismicPolygons, coorIndexes);
             offscreen.startGL();
             offscreen.repaint();
             byte[] polygonData = offscreen.getData();
             offscreen.close();

             seismicCrossplotPolygonPoints = new StsCrossplotPoint[256 * 256];
             byte[][] seismicData = new byte[2][];
             seismicData[0] = seismicVolumes[0].readBytePlaneData(dirNo, dirCoordinate);
             seismicData[1] = seismicVolumes[1].readBytePlaneData(dirNo, dirCoordinate);
             int n = 0;
             int nPoints = 0;
             StsTypeLibrary typeLibrary = crossplot.getTypeLibrary();
             for (int row = 0; row < nRows; row++)
             {
                 for (int col = 0; col < nCols; col++)
                 {
                     int index = StsMath.signedByteToUnsignedInt(polygonData[n]);
                     if (index > 0)
                     {
                         StsColor color = seismicPolygons[index - 1].getStsColor();
                         seismicCrossplotPolygonPoints[nPoints++] = new StsCrossplotPoint(color, row, col, seismicData[0][n], seismicData[1][n]);
                     }
                     n++;
                 }
             }
             seismicCrossplotPolygonPoints = (StsCrossplotPoint[]) StsMath.trimArray(seismicCrossplotPolygonPoints, nPoints);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsCrossplot.processPolygons() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     static public void clearSeismicPoints()
     {
         seismicPolygons = null;
         seismicCrossplotPolygonPoints = null;
         seismicCrossplotPoints = null;
         dirNo = -1;
         dirCoordinate = StsParameters.nullValue;
     }

     private void addSeismicPolygon(StsXPolygon polygon)
     {
         int nPolygons = 0;
         if (seismicPolygons != null) nPolygons = seismicPolygons.length;
         seismicPolygons = (StsXPolygon[]) StsMath.arrayAddElement(seismicPolygons, polygon);
         polygon.setID(nPolygons);
     }

     private void deleteSeismicPolygon(StsXPolygon polygon)
     {
         seismicPolygons = (StsXPolygon[]) StsMath.arrayDeleteElement(seismicPolygons, polygon);
         int nPolygons = seismicPolygons.length;

         // ID used by seismicPolygons is the sequence, so reindex
         for (int n = 0; n < nPolygons; n++)
             seismicPolygons[n].setID(n);

         processSeismicPolygons();
     }

     static public StsXPolygon[] getSeismicPolygons(int dirNo_, float dirCoordinate_)
     {
         if (checkDirAndPlane(dirNo_, dirCoordinate_)) return seismicPolygons;
         return null;
     }

     static private boolean checkDirAndPlane(int dirNo_, float dirCoordinate_)
     {
         return dirNo_ == dirNo && dirCoordinate_ == dirCoordinate;
     }

     static public StsCrossplotPoint[] getSeismicCrossplotPoints(int dirNo_, float dirCoordinate_)
     {
         if (checkDirAndPlane(dirNo_, dirCoordinate_)) return seismicCrossplotPoints;
         return null;
     }

     static public StsCrossplotPoint[] getSeismicCrossplotPolygonPoints(int dirNo_, float dirCoordinate_)
     {
         if (checkDirAndPlane(dirNo_, dirCoordinate_)) return seismicCrossplotPolygonPoints;
         return null;
     }

     private void deletePolygonPoint(StsXPolygon polygon, int index)
     {
         polygon.deletePoint(index);
         //        polygonsChanged();
     }

     /*
         static public void setPolygonColor(Color color, String name)
         {
             if(color == null || name == null)
             {
                 currentColor = null;
                 currentName = null;
             }
             else
             {
                 currentColor = new StsColor(color);
                 currentName = name;
             }
         }
     */
     static public int getPolygonActionMode(StsXPolygon polygon)
     {
         if (polygon == selectedPolygon) return mode;
         else return NONE;
     }

     static public int getPolygonActionType(StsXPolygon polygon)
     {
         if (polygon == selectedPolygon) return type;
         else return NONE;
     }

     static public int getPolygonVertexPickIndex() { return pickIndex; }
 }
