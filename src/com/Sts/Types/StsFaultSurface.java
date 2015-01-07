//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
 import com.Sts.MVC.*;
 import com.Sts.Utilities.*;

 import java.util.*;

public class StsFaultSurface
 {
     StsModelSurface surface;
     StsModel model;
     StsEdgeLoopRadialLinkGrid linkedGrid;
     int nRows, nCols;
     float xInc, yInc;
     int rowMin, rowMax, colMin, colMax;
     StsSectionPoint[][] sectionPointArray; // used only during gapSurfaceGrid for distances to faults

     boolean debug = false;

     // Convenience copies of flags and values
     static public final int MINUS = StsParameters.MINUS;
     static public final int PLUS = StsParameters.PLUS;
     static public final int ROW = StsParameters.ROW;
     static public final int COL = StsParameters.COL;
     static public final int ROWCOL = StsParameters.ROWCOL;
     static public final int NONE = StsParameters.NONE;
     static public final int RIGHT = StsParameters.RIGHT;
     static public final int LEFT = StsParameters.LEFT;
     static public final float largeFloat = StsParameters.largeFloat;

 //    static public final byte NULL_GAP = StsParameters.NULL_GAP; // surface grid is a gap null

     //    static public final byte GAP_NULL = StsParameters.GAP_NULL; // gapPoint value is a null
     static public final byte GAP_GRID = StsParameters.GAP_GRID;    // gridPoint for interpolating gapPoints
     static public final byte GAP_NOT_FILLED = StsParameters.GAP_NOT_FILLED; // gapPoint we need to fill
 //    static public final byte GAP_FILLED = StsParameters.GAP_FILLED; // gapPoint value has been filled
 //    static public final byte GAP_SURF_GRID = StsParameters.GAP_SURF_GRID; // gapPoint value is a surface grid point

     static public final int FAULT = StsParameters.FAULT;
     static public final int AUXILIARY = StsParameters.AUXILIARY;

     static public final float searchDistance = 5.0f;

     public StsFaultSurface(StsModel model, StsModelSurface surface, StsEdgeLoopRadialLinkGrid linkedGrid)
     {
         this.model = model;
         this.surface = surface;
         this.linkedGrid = linkedGrid;

         rowMin = linkedGrid.getRowMin();
         rowMax = linkedGrid.getRowMax();
         colMin = linkedGrid.getColMin();
         colMax = linkedGrid.getColMax();

         nRows = rowMax - rowMin + 1;
         nCols = colMax - colMin + 1;

         xInc = surface.getXInc();
         yInc = surface.getYInc();

         sectionPointArray = new StsSectionPoint[nRows][nCols];
         gapSurfaceGrid();
         if (!Main.isDebug) sectionPointArray = null;
     }

     /** Accessor routines */
     public int getNRows()
     { return nRows; }

     public int getNCols(){ return nCols; }

     private boolean gapSurfaceGrid()
     {
         StsList blockGrids = surface.getBlockGrids();
         int nBlockGrids = blockGrids.getSize();

         for (int n = 0; n < nBlockGrids; n++)
         {
             StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
             if (!gapSurfaceGrid(blockGrid)) return false;
         }
         return true;
     }

     private boolean gapSurfaceGrid(StsBlockGrid blockGrid)
     {
         try
         {
             StsList edges = blockGrid.getEdges();
             int nEdges = edges.getSize();

             for (int n = 0; n < nEdges; n++)
             {
                 StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(n);
                 if (!edge.isFaulted()) continue;

                 StsSection section = edge.getSection();
                 int side = edge.getSide();
                 float faultGap = section.getFaultGap(side);

                 Iterator<StsEdgeLoopRadialGridLink> edgeLinksIterator = blockGrid.getEdgeLoop().getEdgeLinksIterator(edge, PLUS);
                 while (edgeLinksIterator.hasNext())
                 {
                     StsEdgeLoopRadialGridLink link = edgeLinksIterator.next();
                     byte[] insideLinkDirections = link.getLinkedDirections();
                     for (int i = 0; i < insideLinkDirections.length; i++)
                         gapSurfaceGrid(link, section, side, faultGap, insideLinkDirections[i]);
                 }
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsFaultSurface.gapSurfaceGrid(blockGrid) failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     public StsSectionPoint getSectionPoint(int row, int col)
     {
         return sectionPointArray[row - rowMin][col - colMin];
     }

     public void setSectionPoint(int row, int col, StsSectionPoint sectionPoint)
     {
         sectionPointArray[row - rowMin][col - colMin] = sectionPoint;
     }

     private void gapSurfaceGrid(StsEdgeLoopRadialGridLink link, StsSection section, int side, float faultBlockGridSize, int insideLinkIndex)
     {
         int row, col;
         int rowStart, rowEnd, colStart, colEnd;
         int inc;
         float faultGapDistance;
         StsSectionPoint sectionPoint;
         float distance;

         try
         {
             StsEdgeLoopRadialGridLink insideLink = link.getLink(insideLinkIndex);
             StsEdgeLoopRadialGridLink nextLink = insideLink.getLink(insideLinkIndex);
             while (nextLink != null)
             {
                 insideLink = nextLink;
                 nextLink = insideLink.getConnectedLink(insideLinkIndex);
             }

             if (insideLink == null)
             {
                 StsException.systemError("StsFaultSurface.gapSurfaceGrid() failed. " +
                     link.toString() + " has no link inside.");
                 return;
             }

             float rowF = link.getRowF();
             float colF = link.getColF();
             int rowOrCol = StsEdgeLoopRadialGridLink.getLinkIndexRowOrCol(insideLinkIndex);
             if (rowOrCol == NONE) return;
             int direction = StsEdgeLoopRadialGridLink.getLinkIndexDirection(insideLinkIndex);
             if (direction == NONE) return; // direction is PLUS (1) or MINUS (-1) or NONE (0)

             float sectionColF = link.getPoint().getColF(section);

             if (rowOrCol == ROW)
             {
                 row = (int) rowF;
                 faultGapDistance = xInc * faultBlockGridSize;

                 if (direction == PLUS)
                 {
                     colStart = StsMath.ceiling(colF);
                     colEnd = StsMath.floor(insideLink.getColF());
                     if (colEnd < colStart) return;
                     inc = 1;
                 }
                 else
                 {
                     colStart = StsMath.floor(colF);
                     colEnd = StsMath.ceiling(insideLink.getColF());
                     if (colEnd > colStart) return;
                     inc = -1;
                 }

                 col = colStart;
                 int lastColGapped = colStart;
                 while (true)
                 {
                     float gridDistance = Math.abs(colF - col);
                     if (gridDistance <= faultBlockGridSize)
                     {
                         if (Main.debugPoint && row == Main.debugPointRow && col == Main.debugPointCol)
                             StsException.systemDebug("StsFaultSurface.gapSurfaceGrid() called." +
                                 " row: " + row + " col " + col + " norm distance to nearest fault: " + gridDistance * faultGapDistance);
                         // surface.setPointNull(row, col, StsGridPoint.SURF_GAP_SET);
                         lastColGapped = col;
                     }
                     else
                     {
                         sectionPoint = getSetSectionPoint(row, col, section, sectionColF, faultGapDistance, ROW, direction, side);
                         distance = sectionPoint.distance;
                         if(distance < 1.0f) lastColGapped = col;
                         if (distance > searchDistance && sectionPoint.side == side) break;
                         sectionColF = sectionPoint.sectionColF;
                     }
                     if (col == colEnd) break;
                     col += inc;
                 }
                 for(col = colStart; col <= lastColGapped; col++)
                    surface.setPointNull(row, col, StsGridPoint.SURF_GAP_SET);
             }
             else if (rowOrCol == COL)
             {
                 col = (int) colF;
                 faultGapDistance = yInc * faultBlockGridSize;

                 if (direction == PLUS)
                 {
                     rowStart = StsMath.ceiling(rowF);
                     rowEnd = StsMath.floor(insideLink.getRowF());
                     if (rowEnd < rowStart) return;
                     inc = 1;
                 }
                 else
                 {
                     rowStart = StsMath.floor(rowF);
                     rowEnd = StsMath.ceiling(insideLink.getRowF());
                     if (rowEnd > rowStart) return;
                     inc = -1;

                 }

                 row = rowStart;
                 while (true)
                 {
                     float gridDistance = Math.abs(rowF - row);
                     if (gridDistance <= faultBlockGridSize)
                     {
                         if (Main.debugPoint && row == Main.debugPointRow && col == Main.debugPointCol)
                             StsException.systemDebug("StsFaultSurface.gapSurfaceGrid() called." +
                                 " row: " + row + " col " + col + " norm distance to nearest fault: " +
                                 gridDistance * faultGapDistance);
                         surface.setPointNull(row, col, StsGridPoint.SURF_GAP_SET);
                     }
                     else
                     {
                         sectionPoint = getSetSectionPoint(row, col, section, sectionColF, faultGapDistance, ROW, direction, side);
                         distance = sectionPoint.distance;
                         if (distance > searchDistance && sectionPoint.side == side) break;
                         sectionColF = sectionPoint.sectionColF;
                     }
                     if (row == rowEnd) break;
                     row += inc;
                 }
             }
         }
         catch (Exception e)
         {
             StsException.outputException("StsFaultSurface.gapSurfaceGrid() failed.",
                 e, StsException.WARNING);
         }
     }

     /**
      * Get or set sectionPoint at this ROW/COL.  If already set and same section,
      * return the sectionPoint.  If not set or normalized distance to this section
      * is smaller than current at this ROW/COL, set new sectionPoint and return it.
      */
     public StsSectionPoint getSetSectionPoint(int i, int j, StsSection section,
                                               float sectionColF, float faultGapDistance, int rowOrCol, int insideDirection, int correctSide)
     {
         StsSectionPoint currentSectionPoint, sectionPoint;

         /** Check boundary rectangle (if it has been set). */
 //        if(!insideBoundary(i, j)) return null;

         /** Check if point inside grid. */

         if (!insideGrid(i, j))
         {
             StsException.outputException(new StsException(StsException.WARNING,
                 "StsFaultSurface.getSetSectionPoint(i,j)", "i or j out of range:", i, " ", j));
             return null;
 //          return new StsSectionPoint(i, j, section, surface, rowOrCol);
         }

         /** If point already computed for this section, return it. */
         currentSectionPoint = getSectionPoint(i, j);
         if (currentSectionPoint != null && currentSectionPoint.section == section)
             return currentSectionPoint;

         sectionPoint = new StsSectionPoint(i, j, section, surface, rowOrCol);

         /** Find nearest point on section and compute distance to it.
          *  If this fails, treat point as a gapPoint (distance = 0).
          */
         sectionPoint.sectionColF = sectionColF;
         if (!section.computeNearestPoint(sectionPoint))
         {
             sectionPoint.distance = 0.0f;
             return sectionPoint;
         }

         /** Normalize the distance by the faultGapSize */
         if (faultGapDistance > 0.0f) sectionPoint.distance /= faultGapDistance;

         /** If this point is on the wrong side of this section, null the surface point
          *  and don't use it.
          */
 /*
         if(sectionPoint.side != correctSide)
         {
 //            surface.setPointNull(i, j, StsGridPoint.NULL_GAP_SET);
             return sectionPoint;
         }
 */
         /** If normalized distance closer than current one at this point, set it.
          *  Don't use this point if it is off the end of the section.
          */
         if (currentSectionPoint == null || currentSectionPoint.distance > sectionPoint.distance)
             setSectionPoint(i, j, sectionPoint);

         if (Main.debugPoint && i == Main.debugPointRow && j == Main.debugPointCol)
             StsException.systemDebug("StsFaultSurface.getSetSectionPoint() called." +
                 " row: " + i + " col " + j + " norm distance to nearest fault: " + sectionPoint.distance);

         if (sectionPoint.distance <= 1.0f)
 //        if(sectionPoint.distance <= 1.0f && sectionPoint.side == correctSide)
             surface.setPointNull(i, j, StsGridPoint.SURF_GAP_SET);
 //            sectionPoint.setDirectionToSection(insideDirection, correctSide);

         return sectionPoint;
     }

     public boolean insideGrid(int row, int col)
     {
         return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
     }

     public boolean gapBlockGrids()
     {
         StsLinePoints.initialize(Math.max(nRows, nCols));

         StsList blockGrids = surface.getBlockGrids();
         int nBlockGrids = blockGrids.getSize();

         for (int n = 0; n < nBlockGrids; n++)
         {
             StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
             if (!gapBlockGrid(blockGrid)) return false;
         }
         StsLinePoints.delete();
         return true;
     }

     public boolean gapBlockGrid(StsBlockGrid blockGrid)
     {
         try
         {
             StsList edges = blockGrid.getEdges();
             int nEdges = edges.getSize();
             for (int n = 0; n < nEdges; n++)
             {
                 StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(n);
                 StsSection section = edge.getSection();
                 int side = edge.getSide();
                 byte type = section.getType();
                 boolean optimizeOK = section.getFaultGap(side) > 0.0f;

                 Iterator<StsEdgeLoopRadialGridLink> edgeLinksIterator = blockGrid.getEdgeLoop().getEdgeLinksIterator(edge, PLUS);
                 while (edgeLinksIterator.hasNext())
                 {
                     StsEdgeLoopRadialGridLink link = edgeLinksIterator.next();
                     if (!gapBlockGrid(link, section, side, type, blockGrid, optimizeOK)) return false;
                 }
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsFaultSurface.gapBlockGridGrid(blockGrid) failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     private boolean gapBlockGrid(StsEdgeLoopRadialGridLink link, StsSection section,
                                  int side, int type, StsBlockGrid blockGrid, boolean optimizeOK)
     {
         try
         {
             byte[] insideLinkDirections = link.getLinkedDirections();

             int nInsideLinkDirections = insideLinkDirections.length;
             if (nInsideLinkDirections == 0)
             {
                 // if on row and col with no inside connects (corner), set blockGrid type appropriately
                 if (link.isRowAndCol())
                 {
                     int row = link.getRow();
                     int col = link.getCol();

                     if (type == FAULT || surface.isPointNullGap(row, col))
                         blockGrid.setPointType(row, col, GAP_NOT_FILLED);
                     else
                         blockGrid.setPointType(row, col, GAP_GRID);
                 }
             }
             else
             {
                 for (int i = 0; i < nInsideLinkDirections; i++)
                     if (!gapBlockGrid(link, section, side, type, insideLinkDirections[i], blockGrid, optimizeOK))
                         return false;
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsFaultSurface.gapBlockGridGrid(link, edge, blockGrid) failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     private boolean gapBlockGrid(StsEdgeLoopRadialGridLink link, StsSection section, int side,
                                  int type, int insideLinkIndex, StsBlockGrid blockGrid, boolean optimizeOK)
     {
         float rowF, colF;
         int rowOrCol, direction;
         int row, col;
         int rowStart, rowEnd, colStart, colEnd;

         try
         {
             rowF = link.getRowF();
             colF = link.getColF();
             rowOrCol = StsEdgeLoopRadialGridLink.getLinkIndexRowOrCol(insideLinkIndex);
             if (rowOrCol == NONE) return true;
             direction = StsEdgeLoopRadialGridLink.getLinkIndexDirection(insideLinkIndex);
             if (direction == NONE) return true;

             StsEdgeLoopRadialGridLink insideLink = link.getLink(insideLinkIndex);

             if (insideLink == null)
             {
                 StsException.systemError("StsFaultSurface.gapBlockGridGrid() failed. " +
                     link.toString() + " has no link inside.");
                 return true;
             }

             int outsideLinkIndex = StsEdgeLoopRadialGridLink.getOppositeLinkIndex(insideLinkIndex);
             StsEdgeLoopRadialGridLink outsideLink = link.getLink(outsideLinkIndex);
             if (outsideLink != null) outsideLink = outsideLink.getLink(outsideLinkIndex);

             if (rowOrCol == ROW)
             {
                 row = (int) rowF;
                 if (direction == PLUS)
                 {
                     colStart = StsMath.ceiling(colF);
                     colEnd = StsMath.floor(insideLink.getColF());
                     StsLinePoints.initializeLine(ROW, row, optimizeOK);
                     for (col = colStart; col <= colEnd; col++)
                         if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                     StsLinePoints.setGoodBlockGridPoints(blockGrid);

                     if (outsideLink != null)
                     {
                         colStart--;
                         colEnd = StsMath.ceiling(outsideLink.getColF());
                         colEnd = Math.max(colEnd, blockGrid.getPointColMin());
                         /*
                             if(type == FAULT)
                             {
                                 for(col = colStart; col >= colEnd; col--)  // clear outside points
                                     blockGrid.setPointTypeIfNot(row, col, GAP_NULL, GAP_GRID);
                             }
                             else if(type == AUXILIARY)
                         */
                         if (type == AUXILIARY)
                         {
                             StsLinePoints.initializeLine(ROW, row, false);
                             for (col = colStart; col >= colEnd; col--)  // clear outside points
                                 if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                             StsLinePoints.setGoodBlockGridPoints(blockGrid);
 //                            if(colStart >= colEnd)
 //                                blockGrid.setPointTypeIf(row, colStart, GAP_GRID, GAP_SURF_GRID);
                         }
                     }
                 }
                 else
                 {
                     colStart = StsMath.floor(colF);
                     colEnd = StsMath.ceiling(insideLink.getColF());
                     StsLinePoints.initializeLine(ROW, row, optimizeOK);
                     for (col = colStart; col >= colEnd; col--)
                         if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                     StsLinePoints.setGoodBlockGridPoints(blockGrid);

                     if (outsideLink != null)
                     {
                         colStart++;
                         colEnd = StsMath.floor(outsideLink.getColF());
                         colEnd = Math.min(colEnd, blockGrid.getPointColMax());
                         /*
                             if(type == FAULT)
                             {
                                 for(col = colStart; col <= colEnd; col++)  // clear outside points
                                     blockGrid.setPointTypeIfNot(row, col, GAP_NULL, GAP_GRID);
                             }
     //                        else if(type == AUXILIARY)
                         */
                         if (type == AUXILIARY)
                         {
                             StsLinePoints.initializeLine(ROW, row, false);
                             for (col = colStart; col <= colEnd; col++)  // clear outside points
                                 if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                             StsLinePoints.setGoodBlockGridPoints(blockGrid);
 //                            if(colStart <= colEnd)
 //                                blockGrid.setPointTypeIf(row, colStart, GAP_GRID, GAP_SURF_GRID);
                         }
                     }
                 }
             }
             else if (rowOrCol == COL)
             {
                 col = (int) colF;
                 if (direction == PLUS)
                 {
                     rowStart = StsMath.ceiling(rowF);
                     rowEnd = StsMath.floor(insideLink.getRowF());
                     StsLinePoints.initializeLine(COL, col, optimizeOK);
                     for (row = rowStart; row <= rowEnd; row++)
                         if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                     StsLinePoints.setGoodBlockGridPoints(blockGrid);

                     if (outsideLink != null)
                     {
                         rowStart--;
                         rowEnd = StsMath.ceiling(outsideLink.getRowF());
                         rowEnd = Math.max(rowEnd, blockGrid.getPointRowMin());
                         /*
                             if(type == FAULT)
                             {
                                 for(row = rowStart; row >= rowEnd; row--)  // clear outside points
                                     blockGrid.setPointTypeIfNot(row, col, GAP_NULL, GAP_GRID);
                             }
     //                        else if(type == AUXILIARY)
                         */
                         if (type == AUXILIARY)
                         {
                             StsLinePoints.initializeLine(COL, col, false);
                             for (row = rowStart; row >= rowEnd; row--)  // clear outside points
                                 if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                             StsLinePoints.setGoodBlockGridPoints(blockGrid);
 //                            if(rowStart >= rowEnd)
 //                                blockGrid.setPointTypeIf(rowStart, col, GAP_GRID, GAP_SURF_GRID);
                         }
                     }
                 }
                 else
                 {
                     rowStart = StsMath.floor(rowF);
                     rowEnd = StsMath.ceiling(insideLink.getRowF());
                     StsLinePoints.initializeLine(COL, col, optimizeOK);
                     for (row = rowStart; row >= rowEnd; row--)
                         if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                     StsLinePoints.setGoodBlockGridPoints(blockGrid);

                     if (outsideLink != null)
                     {
                         rowStart++;
                         rowEnd = StsMath.floor(outsideLink.getRowF());
                         rowEnd = Math.min(rowEnd, blockGrid.getPointRowMax());
                         /*
                             if(type == FAULT)
                             {
                                 for(row = rowStart; row <= rowEnd; row++)  // clear outside points
                                     blockGrid.setPointTypeIfNot(row, col, GAP_NULL, GAP_GRID);
                             }
     //                        else if(type == AUXILIARY)
                         */
                         if (type == AUXILIARY)
                         {
                             StsLinePoints.initializeLine(COL, col, false);
                             for (row = rowStart; row <= rowEnd; row++)  // clear outside points
                                 if (!setGridOrFill(row, col, blockGrid, section, side)) break;
                             StsLinePoints.setGoodBlockGridPoints(blockGrid);
 //                            if(rowStart <= rowEnd)
 //                                blockGrid.setPointTypeIf(rowStart, col, GAP_GRID, GAP_SURF_GRID);
                         }
                     }
                 }
             }

             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsFaultSurface.gapBlockGridGrid() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     /**
      * For this point, flag as GAP_NOT_FILLED (gapPoint) or GAP_GRID (use in
      * interpolation).  Return true if we can continue searching.  Process
      * stops when we have found optimum fit for a set of GAP_GRID points we have
      * have found.
      */
     private boolean setGridOrFill(int row, int col, StsBlockGrid blockGrid, StsSection section, int side)
     {
         /** If already GAP_NOT_FILLED we have determined its in the gap
          *  or a poor-quality point, so return and check next one.
          */
         if (blockGrid.getPointType(row, col) == GAP_NOT_FILLED) return true;

         /** If surface has null gap, set point to GAP_NOT_FILLED so it
          *  will be subsequently computed and filled.
          */
         if (surface.isPointNullGap(row, col))
         {
             blockGrid.setPointType(row, col, GAP_NOT_FILLED);
             return true;
         }

         /** This is a goodPoint candidate.  Add it to linePoints and see if it
          *  improves current fit; return true if it does, false if not.
          */

         float z = surface.getZorT(row, col);
         return StsLinePoints.addPointBetterFit(row, col, z);
     }
 }
