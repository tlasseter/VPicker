//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
 import com.Sts.Interfaces.*;
 import com.Sts.MVC.*;
 import com.Sts.Utilities.*;
 import com.Sts.Utilities.Geometry.*;

 import javax.media.opengl.*;
 import javax.media.opengl.glu.*;
 import java.awt.event.*;
 import java.util.*;

public class StsPolygon implements IPolygon, Comparable<StsPolygon>
 {
     public float[][] pntNrmls;
     public float[] fValues;
     public StsGridRowCol[] gridRowCols;
     public int nPoints;
     public int row = -1, col = -1, layer = -1, nBlock = -1;
     /**
      * type can be: SURFACE (horizon polygon), SECTION (section polygon), GRID (on section bounded by grid column),
      * LAYER_RIGHT and LAYER_LEFT (section subdivision of grid polygon into right and left side layer subpolygons).
      */
     public byte type;
     StsList links;
     private StsList subPolygons;
     public StsColor stsColor = null;
     public float[] center;
     // public float area;

     static StsConcavePolygon concavePolygon = null; // extends GLUtesselatorCallbackAdapter
     static GLUtessellator tesselator = null;
     static boolean drawConvexOnly = false;
     static boolean drawingConcave = true;

     static boolean debugPolygons = true;
     static int debugPolygonRow = 9;
     static int debugPolygonCol = 97;
     static int debugPolygonMinLayer = 8;
     static int debugPolygonMaxLayer = 9;
     static boolean debug = false;

     static final float largeFloat = StsParameters.largeFloat;
     static public final int MAX_POINTS = 1000;
     static public final int EXCESSIVE_POINTS = 100;

     /**
      * polygons are organized in a hierrarchical structure; trim loop on a section creates SECTION polygons
      * which are cut by grid lines forming GRID subpolygons, which are cut by layers on each side of section
      * forming LAYER_LEFT and LAYER_RIGHT polygons.
      */
     static public final byte SURFACE = 1;
     static public final byte SECTION = 2;
     static public final byte GRID = 3;
     static public final byte LAYER_RIGHT = 4;
     static public final byte LAYER_LEFT = 5;

     static public final String[] typeStrings = new String[]{"NONE", "SURFACE", "SECTION", "GRID", "LAYER_RIGHT", "LAYER_LEFT"};

     static final int ROW = StsParameters.ROW;
     static final int COL = StsParameters.COL;
     static final int ROWCOL = StsParameters.ROWCOL;

     static final int PLUS = StsParameters.PLUS;
     static final int MINUS = StsParameters.MINUS;
     static final int NONE = StsParameters.NONE;
     static final float roundOff = StsParameters.roundOff;

     public StsPolygon()
     {
     }

     public StsPolygon(StsEdgeLoopRadialGridLink firstLink, byte type)
     {
         addLink(firstLink);
         this.type = type;
         initializeDebug();
     }

     public StsPolygon(StsList links, byte type, StsSurfaceGridable grid)
     {
         this.type = type;

         StsEdgeLoopRadialGridLink link;
         if (links == null) return;
         int nLinks = links.getSize();
         for (int n = 0; n < nLinks; n++)
         {
             link = (StsEdgeLoopRadialGridLink) links.getElement(n);
             addLink(link);
         }
         StsEdgeLoopRadialGridLink link0 = (StsEdgeLoopRadialGridLink) links.getFirst();
         StsEdgeLoopRadialGridLink link1 = (StsEdgeLoopRadialGridLink) links.getSecond();
         int[] rowCol = link0.getLowerLeftRowCol(link1, type, grid);
         if (rowCol == null) return;
         this.row = rowCol[0];
         this.col = rowCol[1];
         initializeDebug();
     }

     public StsPolygon(int row, int col)
     {
         this.row = row;
         this.col = col;
         initializeDebug();
     }

     public int getRow(){ return row; }

     public int getCol(){ return col; }

     public void setType(byte type){ this.type = type; }

     public int getType(){ return type; }

     public void setPntNrmls(float[][] pntNrmls)
     {
         this.pntNrmls = pntNrmls;
         if (pntNrmls != null) nPoints = pntNrmls.length;
     }

     public void initializeLayer(int layer, byte layerSide, StsZone zone, float[][] points, int blockNumber)
     {
         layer = zone.checkLimitLayerNumber(layer);
         this.layer = layer;
         type = layerSide;
         // stsColor = zone.getLayerColor(layer);
         setPntNrmls(points);
         nBlock = blockNumber;
         if (debug && layer >= debugPolygonMinLayer && layer <= debugPolygonMaxLayer)
             StsException.systemDebug(this, "setLayer", "polygon: " + toString());
     }

     public void setLayer(float layerF, byte layerSide, StsZone zone)
     {
         layer = zone.checkLimitLayerNumber(layerF);
         type = layerSide;
         stsColor = zone.getLayerColor(layer);
     }

     public void setBlockNumber(int blockNumber)
     {
         this.nBlock = blockNumber;
     }

     public void setStsColor(StsColor color){ stsColor = color; }

     public StsList getLinks(){ return links; }

     static public void toggleDrawConvex(ItemEvent e)
     {
         drawConvexOnly = (e.getStateChange() == ItemEvent.SELECTED);
     }

     static public boolean getDrawConvexOnly(){ return drawConvexOnly; }

     private void initializeDebug()
     {
         debug = debugPolygons && debugPolygonRow == row && debugPolygonCol == col;
         if (debug)
         {
             StsException.systemDebug("StsPolygon initializeDebug called. " +
                 "at row: " + row + " col: " + col);
         }
     }

     public StsEdgeLoopRadialGridLink getLinkAfter(StsEdgeLoopRadialGridLink link)
     {
         int index = links.getIndex(link);
         if (index == StsList.NO_MATCH) return null;
         return (StsEdgeLoopRadialGridLink) links.getElement((index + 1) % links.getSize());
     }

     public void addLink(StsEdgeLoopRadialGridLink link)
     {
         if (links == null) links = new StsList(4, 4);
         links.add(link);
     }

     public void addLinks(StsEdgeLoopRadialGridLink[] addLinks, boolean forwardOrder)
     {
         int n;

         int nLinks = addLinks.length;
         if (links == null) links = new StsList(nLinks);
         if (forwardOrder)
             for (n = 0; n < nLinks; n++)
                 links.add(addLinks[n]);
         else
             for (n = nLinks - 1; n >= 0; n--)
                 links.add(addLinks[n]);
     }

     public boolean checkOK(StsSurfaceGridable grid)
     {
         if (links == null) return false;
         int nLinks = links.getSize();
         if (nLinks < 3) return false;
         /*
                 for(int n = 0; n < nLinks; n++)
                 {
                     StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getElement(n);
                     if(!link.usageOK()) return false;
                 }
         */
         return true;
     }

     public void setLowerLeftRowCol(StsSurfaceGridable grid)
     {
         int[] rowCol = getLowerLeftRowCol(grid);
         row = rowCol[0];
         col = rowCol[1];
     }

     public void setLowerLeftRowCol(StsSurfaceGridable grid, byte polygonType)
     {
         int[] rowCol = getLowerLeftRowCol(grid, polygonType);
         row = rowCol[0];
         col = rowCol[1];
     }

     public int[] getLowerLeftRowCol(StsSurfaceGridable grid)
     {
         if (grid instanceof StsSection)
         {
             int[] rowCol = getMinimumRowCol(grid, StsPolygon.SECTION);
             adjustLowerLeftRowCol(rowCol, (StsSection) grid);
             return rowCol;
         }
         else
             return getMinimumRowCol(grid, StsPolygon.GRID);
     }

     public int[] getLowerLeftRowCol(StsSurfaceGridable grid, byte polygonType)
     {
         if (grid instanceof StsSection)
         {
             int[] rowCol = getMinimumRowCol(grid, polygonType);
             adjustLowerLeftRowCol(rowCol, (StsSection) grid);
             return rowCol;
         }
         else
             return getMinimumRowCol(grid, StsPolygon.GRID);
     }

     public int[] getMinimumRowCol(StsSurfaceGridable grid, byte polygonType)
     {
         int row = StsParameters.largeInt;
         int col = StsParameters.largeInt;
         int nLinks = links.getSize();
         for (int n = 0; n < nLinks; n++)
         {
             StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) links.getElement(n);
             StsRowCol rowCol = link.point.getRowCol(grid, polygonType);
             row = Math.min(row, StsMath.floor(rowCol.getRowF()));
             col = Math.min(col, StsMath.floor(rowCol.getColF()));
         }
         return new int[]{row, col};
     }

     public void adjustLowerLeftRowCol(int[] rowCol, StsSection section)
     {
         byte geometryType = section.getGeometry().getGeometryType();
         if (geometryType == StsSection.GEOM_ROW_PLUS)
             rowCol[0]--;
         else if (geometryType == StsSection.GEOM_COL_MINUS)
             rowCol[1]--;
     }

     public void computePoints(StsSurfaceGridable grid, boolean isPlanar)
     {
         int nLinks = links.getSize();
         nPoints = nLinks;
         pntNrmls = new float[nLinks][];
         gridRowCols = new StsGridRowCol[nLinks];
         if (isPlanar)
         {
             for (int n = 0; n < nLinks; n++)
             {
                 StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) links.getElement(n);
                 pntNrmls[n] = link.getXYZ();
                 gridRowCols[n] = link.getPoint().getGridRowCol();
             }
         }
         else
         {
             for (int n = 0; n < nLinks; n++)
             {
                 StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) links.getElement(n);
                 pntNrmls[n] = link.getVertexAndNormal(grid);
                 gridRowCols[n] = link.getPoint().getGridRowCol();
             }
         }
         if (!Main.isDebug) links = null;
     }

     // starting from the first point, get the last point which is on the same edge
     public StsEdgeLoopRadialGridLink getNextEdgeLink(StsEdgeLoopRadialGridLink firstLink)
     {
         StsEdgeLoopRadialGridLink link, nextLink;

         if (links == null) return firstLink.getNextEdgeLink();
         int nLinks = links.getSize();
         if (nLinks == 0) return firstLink.getNextEdgeLink();

         int firstIndex = links.getIndex(firstLink);
         if (firstIndex == StsList.NO_MATCH) return firstLink.getNextEdgeLink();
         nextLink = (StsEdgeLoopRadialGridLink) links.getElement(firstIndex);
         for (int n = firstIndex + 1; n < nLinks; n++)
         {
             link = nextLink;
             nextLink = (StsEdgeLoopRadialGridLink) links.getElement(n);
             if (nextLink != link.getNextEdgeLink()) return link;
         }
         return firstLink.getNextEdgeLink();
     }

     public void setSubPolygons(StsList subPolygons)
     {
         this.subPolygons = subPolygons;
     }

     public void checkSize()
     {
         int nLinks = links.getSize();
         if (nLinks < EXCESSIVE_POINTS) return;
         StsException.systemError("Number of links in polygon: " + getLabel() +
             " is excessive: " + nLinks);
     }

     // Points form a clockwise polygon.  Find topZ and botZ point indices.
     // Compute the f values for each point: f ranges from 0 to the number of subZones.
     // If two points are even at top, these form the top points of the ladder.
     // If two points are even at the bot, these form the bot points of the ladder.
     // If either end is not even, then the extreme point is used on both sides of the ladder.
     // For the ladder, insert rungs at integral values of f.
     // Make subPolygons between each rung and above the top rung and below the bottom rung.
     public void constructZoneSideLayerPolygons(StsZoneBlock zoneBlock, byte layerSide)
     {
         try
         {
             // If GRID, subPolygons are subZones or layers: delete because
             // we are recomputing.
             if (type == GRID && subPolygons != null) subPolygons = null;
             int blockNumber = zoneBlock.getBlock().getIndex();
             // if subPolygons aren't null, this is a SECTION polygon, i.e., it has GRID subPolygons so divide these GRID subPolygons into layerCells
             if (subPolygons != null)
             {
                 int nSubPolygons = subPolygons.getSize();
                 for (int n = 0; n < nSubPolygons; n++)
                 {
                     StsPolygon subPolygon = (StsPolygon) subPolygons.getElement(n);
                     if (subPolygon != null) subPolygon.constructZoneSideLayerPolygons(zoneBlock, layerSide);
                 }
                 nBlock = blockNumber;
                 return;
             }

             if (pntNrmls == null) return;

             // rowCol values are the float grid row & col numbers for each polygon point
             checkComputeRowCols(zoneBlock);

             // f values are scaled from 0 at top horizon to nSubZones at bot horizon
             checkComputeLayerFValues(zoneBlock);

             // construct temp class instances with pntNrml and fValue for each point
             PointPlusF[] points = new PointPlusF[nPoints];
             float fMin = StsParameters.largeFloat;
             float fMax = -StsParameters.largeFloat;
             for (int n = 0; n < nPoints; n++)
             {
                 float f = fValues[n];
                 points[n] = new PointPlusF(pntNrmls[n], f, n);
                 fMin = Math.min(f, fMin);
                 fMax = Math.max(f, fMax);
             }
             // insert integer f points in between point pairs
             PointPlusF p0, p1;
             int if0, if1;
             int n0, n1;
             TreeSet leftRungSet = new TreeSet();
             TreeSet riteRungSet = new TreeSet();

             p1 = points[0];
             for (n0 = 0; n0 < nPoints; n0++)
             {
                 p0 = p1;
                 n1 = (n0 + 1) % nPoints;
                 p1 = points[n1];

                 if (p0.f < p1.f)
                 {
                     if0 = StsMath.ceiling(p0.f);
                     if1 = StsMath.floor(p1.f);
                     for (int i = if0; i <= if1; i++)
                         riteRungSet.add(interpolatePointPlusF(p0, p1, i, n0));
                 }
                 else if (p0.f > p1.f)
                 {
                     if0 = StsMath.floor(p0.f);
                     if1 = StsMath.ceiling(p1.f);
                     for (int i = if0; i >= if1; i--)
                         leftRungSet.add(interpolatePointPlusF(p0, p1, i, n0));
                 }
                 /*
                 else // point heights are even: add first to left and second to right
                 {
                     int rung = (int)p0.f;
                     p0.rung = rung;
                     p1.rung = rung;
                     leftRungSet.add(p0);
                     riteRungSet.add(p1);
                 }
                 */
             }

             int nLeftRungs = leftRungSet.size();
             int nRiteRungs = riteRungSet.size();
             if (nLeftRungs != nRiteRungs)
             {
                 StsException.systemError("constructZoneSidePolygons() failed." +
                     " Number of rung sides not agree. Left: " + nLeftRungs + " Rite: " + nRiteRungs);
                 return;
             }

             int nRungs = nLeftRungs;

             PointPlusF topLeftRung, topRiteRung, botLeftRung, botRiteRung;

             // eliminate top or bot rung if both sides are same point
             if (nRungs > 0)
             {
                 topLeftRung = (PointPlusF) leftRungSet.first();
                 topRiteRung = (PointPlusF) riteRungSet.first();
                 if (topLeftRung.indexF == topRiteRung.indexF)
                 {
                     leftRungSet.remove(topLeftRung);
                     riteRungSet.remove(topRiteRung);
                     nRungs--;

                 }
             }
             if (nRungs > 0)
             {
                 botLeftRung = (PointPlusF) leftRungSet.last();
                 botRiteRung = (PointPlusF) riteRungSet.last();
                 if (botLeftRung.indexF == botRiteRung.indexF)
                 {
                     leftRungSet.remove(botLeftRung);
                     riteRungSet.remove(botRiteRung);
                     nRungs--;
                 }
             }

             StsZone zone = zoneBlock.getZone();
             if (nRungs == 0)
             {
                 setLayer(fMin, layerSide, zone);
                 return;
             }

             Object[] leftRungs = leftRungSet.toArray();
             Object[] riteRungs = riteRungSet.toArray();

             StsPolygon polygon;
             int nSubPoints, nSubPoint;
             float[][] subPoints;

             boolean topEven = ((PointPlusF) leftRungs[0]).f == fMin;
             boolean botEven = ((PointPlusF) leftRungs[nRungs - 1]).f == fMax;

             int nPolygons = nRungs + 1;
             if (topEven) nPolygons--;
             if (botEven) nPolygons--;
             subPolygons = new StsList(nPolygons, 2);

             botLeftRung = (PointPlusF) leftRungs[0];
             botRiteRung = (PointPlusF) riteRungs[0];
             int[] gridRowCol = StsRowCol.getLowerLeftRowCol(gridRowCols);
             int gridRow = gridRowCol[0];
             int gridCol = gridRowCol[1];
             // build polygon above topRung
             if (!topEven)
             {
                 nSubPoints = 2 + numberIntermediatePoints(botLeftRung, botRiteRung);
                 if (nSubPoints > 2)
                 {
                     subPoints = new float[nSubPoints][];

                     nSubPoint = 0;
                     subPoints[nSubPoint++] = botLeftRung.pntNrml;
                     nSubPoint = copyPntNrmls(pntNrmls, botLeftRung, botRiteRung, subPoints, nSubPoint);
                     subPoints[nSubPoint++] = botRiteRung.pntNrml;
                     polygon = new StsPolygon(gridRow, gridCol);
                     int layer = botLeftRung.rung - 1;
                     polygon.initializeLayer(layer, layerSide, zone, subPoints, blockNumber);
                     subPolygons.add(polygon);
                 }
             }

             // build polygons between rungs
             for (int n = 1; n < nRungs; n++)
             {
                 topLeftRung = botLeftRung;
                 topRiteRung = botRiteRung;
                 botLeftRung = (PointPlusF) leftRungs[n];
                 botRiteRung = (PointPlusF) riteRungs[n];

                 nSubPoints = 4 + numberIntermediatePoints(botLeftRung, topLeftRung) +
                     numberIntermediatePoints(topRiteRung, botRiteRung);

                 subPoints = new float[nSubPoints][];

                 nSubPoint = 0;
                 subPoints[nSubPoint++] = botLeftRung.pntNrml;
                 nSubPoint = copyPntNrmls(pntNrmls, botLeftRung, topLeftRung, subPoints, nSubPoint);
                 subPoints[nSubPoint++] = topLeftRung.pntNrml;
                 subPoints[nSubPoint++] = topRiteRung.pntNrml;
                 nSubPoint = copyPntNrmls(pntNrmls, topRiteRung, botRiteRung, subPoints, nSubPoint);
                 subPoints[nSubPoint++] = botRiteRung.pntNrml;

                 polygon = new StsPolygon(gridRow, gridCol);
                 int layer = topLeftRung.rung;
                 polygon.initializeLayer(layer, layerSide, zone, subPoints, blockNumber);
                 subPolygons.add(polygon);
             }

             //build polygon below bottom rung
             if (!botEven)
             {
                 topLeftRung = botLeftRung;
                 topRiteRung = botRiteRung;

                 nSubPoints = 2 + numberIntermediatePoints(topRiteRung, topLeftRung);
                 if (nSubPoints > 2)
                 {
                     subPoints = new float[nSubPoints][];

                     nSubPoint = 0;
                     subPoints[nSubPoint++] = topLeftRung.pntNrml;
                     subPoints[nSubPoint++] = topRiteRung.pntNrml;
                     nSubPoint = copyPntNrmls(pntNrmls, topRiteRung, topLeftRung, subPoints, nSubPoint);

                     polygon = new StsPolygon(gridRow, gridCol);
                     int layer = topLeftRung.rung;
                     polygon.initializeLayer(layer, layerSide, zone, subPoints, blockNumber);
                     subPolygons.add(polygon);
                 }
             }
             type = GRID;
             this.nBlock = blockNumber;
         }
         catch (Exception e)
         {
             StsException.outputException("StsPolygon.constructZoneSideSubPolygons() failed.",
                 e, StsException.WARNING);
         }
     }

     static public int[] getPolygonRowCol(float[] firstRowCol, float[] nextRowCol)
     {
         float rowF = firstRowCol[0];
         int row = (int) rowF;
         float colF = firstRowCol[1];
         int col = (int) colF;
         boolean isRow = (row == rowF);
         boolean isCol = (col == colF);
         if (!isRow && !isCol) return new int[]{row, col};
         float nextRowF = nextRowCol[0];
         float nextColF = nextRowCol[1];
         if (isRow)
         {
             if (nextRowF < rowF) row--;
             else if (nextRowF == rowF && nextColF > colF) row--;
         }
         if (isCol)
         {
             if (nextColF < colF) col--;
             else if (nextColF == colF && nextRowF < rowF) col--;
         }
         return new int[]{row, col};
     }

     private void checkComputeRowCols(StsZoneBlock zoneBlock)
     {
         if (gridRowCols != null) return;
         StsBlockGrid grid = zoneBlock.getTopGrid();
         gridRowCols = new StsGridRowCol[nPoints];
         for (int n = 0; n < nPoints; n++)
             gridRowCols[n] = new StsGridRowCol(grid, pntNrmls[n]);
     }

     private void checkComputeLayerFValues(StsZoneBlock zoneBlock)
     {
         // f values are scaled from 0 at top horizon to nLayers at bot horizon
         if (fValues != null) return;
         fValues = new float[nPoints];
         for (int n = 0; n < nPoints; n++)
             fValues[n] = zoneBlock.getLayerF(gridRowCols[n].rowF, gridRowCols[n].colF, pntNrmls[n][2]);
     }

     public int[] getLayerRange(StsZoneBlock zoneBlock)
     {
         checkComputeRowCols(zoneBlock);
         checkComputeLayerFValues(zoneBlock);
         float fMin = fValues[0];
         float fMax = fValues[0];
         for (int n = 1; n < nPoints; n++)
         {
             float f = fValues[n];
             if (f < fMin)
                 fMin = f;
             else if (f > fMax)
                 f = fMax;
         }
         int layerMin = StsMath.floor(fMin);
         int layerMax = StsMath.below(fMax);
         return new int[]{layerMin, layerMax};
     }

     private int copyPntNrmls(float[][] pntNrmls, PointPlusF p0, PointPlusF p1, float[][] subPoints, int nSubPoint)
     {
         try
         {
             float fi0 = p0.indexF;
             float fi1 = p1.indexF;
             if (fi0 > fi1) fi1 += nPoints;

             int i0 = StsMath.above(fi0);
             int i1 = StsMath.below(fi1);

             for (int i = i0; i <= i1; i++)
                 subPoints[nSubPoint++] = pntNrmls[i % nPoints];
             return nSubPoint;
         }
         catch (Exception e)
         {
             StsException.outputException("StsPolygon.copyPntNrmls() failed.",
                 e, StsException.WARNING);
             return nSubPoint;
         }
     }

     private int numberIntermediatePoints(PointPlusF p0, PointPlusF p1)
     {
         float fi0 = p0.indexF;
         float fi1 = p1.indexF;
         if (fi0 > fi1) fi1 += nPoints;
         return Math.max(0, StsMath.below(fi1) - StsMath.above(fi0) + 1);
     }

     private PointPlusF interpolatePointPlusF(PointPlusF p0, PointPlusF p1, int rung, int n0)
     {
         float df = p1.f - p0.f;
         float ff = (rung - p0.f) / df;
         if (ff == 0.0f)
         {
             p0.rung = rung;
             return p0;
         }
         else if (ff == 1.0f)
         {
             p1.rung = rung;
             return p1;
         }

         float[] pntNrml = StsMath.interpolate(p0.pntNrml, p1.pntNrml, ff);
         float indexF = n0 + ff;
         float f = p0.f + ff * df;

         PointPlusF pointF = new PointPlusF(pntNrml, f, indexF);
         pointF.rung = rung;
         return pointF;
     }

     // if there are no points between these or they are at same f, then they are "even"
     private boolean isEven(PointPlusF p0, PointPlusF p1)
     {
         if (p0.rung != p1.rung) return false;

         float fi0 = p0.indexF;
         float fi1 = p1.indexF;
         if (fi0 > fi1) fi1 += nPoints;

         int i0 = StsMath.above(fi0);
         int i1 = StsMath.below(fi1);

         for (int i = i0; i <= i1; i++)
             if (fValues[i % nPoints] != p0.f) return false;

         return true;

         //        float dIndexF = Math.abs(p1.indexF - p0.indexF);
         //        return dIndexF == 1.0f || dIndexF == nPoints - 1.0f;
     }

     public StsList getSubPolygons()
     {
         return subPolygons;
     }

     public boolean hasSubPolygons()
     {
         return subPolygons != null && subPolygons.getSize() > 0;
     }

     public Iterator<StsPolygon> getSubPolygonIterator()
     {
         ArrayList<StsPolygon> polygons = new ArrayList<StsPolygon>();
         int nPolygons = subPolygons.getSize();
         for (int n = 0; n < nPolygons; n++)
             polygons.add((StsPolygon) subPolygons.getElement(n));
         return polygons.iterator();
     }

     class PointPlusF implements Comparable
     {
         public float[] pntNrml; // float[6]: point xyz and normal xyz
         public float f; // fractional distance from top to bot horizons (between 0 and 1)
         public float indexF; // position in sequence of polygon points
         public int rung; // rung number: corresponds to subZone number

         PointPlusF(float[] pntNrml, float f, int n)
         {
             this(pntNrml, f, (float) n);
         }

         PointPlusF(float[] pntNrml, float f, float indexF)
         {
             this.pntNrml = pntNrml;
             this.f = f;
             this.indexF = indexF;
             rung = -1;
         }

         public int compareTo(Object other)
         {
             return rung - ((PointPlusF) other).rung;
         }

         public int getPrevIndex()
         {
             if (indexF == 0.0f) return nPoints - 1;
             else return StsMath.below(indexF);
         }

         public int getNextIndex()
         {
             if (indexF == (float) nPoints - 1) return 0;
             else return StsMath.above(indexF);
         }
     }

     public boolean draw(GL gl, GLU glu, boolean isPlanar, float[] normal, boolean displayLayers, boolean displayProperties,
                         StsPropertyVolume propertyVolume, boolean propertyChanged)
     {
         if (normal == null)
             normal = getFirstPointNormal();
         if (subPolygons != null && displayLayers)
         {
             int nSubPolygons = subPolygons.getSize();
             for (int n = 0; n < nSubPolygons; n++)
             {
                 StsPolygon polygon = (StsPolygon) subPolygons.getElement(n);
                 if (polygon != null)
                     if (!polygon.draw(gl, glu, isPlanar, normal, true, displayProperties, propertyVolume, propertyChanged))
                         polygon = null;
             }
             return true;
         }
         else
         {
             if (pntNrmls == null) return false;

             if (displayProperties)
             {
                 if (propertyChanged && propertyVolume != null)
                 {
                     stsColor = setPropertyColor(propertyVolume);
                     if (stsColor == null) stsColor = getLayerColor();
                 }
                 if (stsColor != null) stsColor.setGLColor(gl);
             }
             else if (displayLayers)
             {
                 stsColor = getLayerColor();
                 stsColor.setGLColor(gl);
             }

             if (nPoints <= 3 || drawConvexOnly)
             {
                 drawConvex(gl, isPlanar, normal);
                 return true;
             }
             else
                 return drawConcave(gl, glu, isPlanar, normal);
             //                if(!ok) drawConvex(gl, isPlanar, normal);
         }
     }

     private StsColor getLayerColor()
     {
         int colorIndex = 0;
         if (layer > 0)
             colorIndex = layer % 32;
         return StsColor.colors32[colorIndex];
     }

     public boolean drawContours(GL gl, GLU glu, boolean isPlanar, float[] normal, boolean drawSubPolygons)
     {
         if (normal == null)
             normal = getFirstPointNormal();
         if (subPolygons != null && drawSubPolygons)
         {
             int nSubPolygons = subPolygons.getSize();
             for (int n = 0; n < nSubPolygons; n++)
             {
                 StsPolygon polygon = (StsPolygon) subPolygons.getElement(n);
                 if (polygon != null)
                     if (!polygon.drawContours(gl, glu, isPlanar, normal, true)) polygon = null;
             }
             return true;
         }
         else
         {
             if (pntNrmls == null) return false;

             if (stsColor != null) stsColor.setGLColor(gl);

             if (nPoints <= 3 || drawConvexOnly)
             {
                 drawConvex(gl, isPlanar, normal);
                 return true;
             }
             else
                 return drawConcave(gl, glu, isPlanar, normal);
             //                if(!ok) drawConvex(gl, isPlanar, normal);
         }
     }

     private float[] getFirstPointNormal()
     {
         float[] normal = new float[3];
         if (pntNrmls[0].length == 3) return normal;
         for (int n = 0; n < 3; n++)
             normal[n] = pntNrmls[0][n + 3];
         return normal;
     }

     private void drawConvex(GL gl, boolean isPlanar, float[] normal)
     {
         gl.glBegin(GL.GL_POLYGON);

         if (isPlanar)
         {
             gl.glNormal3fv(normal, 0);
             for (int n = 0; n < nPoints; n++)
             {
                 gl.glNormal3fv(normal, 0);
                 gl.glVertex3fv(pntNrmls[n], 0);
             }
         }
         else
         {
             for (int n = 0; n < nPoints; n++)
             {
                 float[] pn = pntNrmls[n];
                 gl.glNormal3f(pn[3], pn[4], pn[5]);
                 gl.glVertex3fv(pn, 0);
             }
         }

         gl.glEnd();
     }

     private boolean drawConcave(GL gl, GLU glu, boolean isPlanar, float[] normal)
     {
         drawingConcave = true;
         if (tesselator == null)
         {
             tesselator = glu.gluNewTess();
             concavePolygon = new StsConcavePolygon(gl, glu);
         }
         String name = toString();
         concavePolygon.initialize(tesselator, this, debug, isPlanar, normal, name);
         glu.gluBeginPolygon(tesselator);
         int nCoords = 6;
         if (isPlanar) nCoords = 3;
         for (int n = 0; n < nPoints; n++)
         {
             double[] coords = new double[nCoords];
             for (int i = 0; i < nCoords; i++)
                 coords[i] = (double) pntNrmls[n][i];
             glu.gluTessVertex(tesselator, coords, 0, coords);
         }
         glu.gluEndPolygon(tesselator);
         return drawingConcave;
     }

     public void drawConcaveFailed(int error)
     {
         drawingConcave = false;
     }

     public void drawConcaveFailed(String errorString)
     {
         drawingConcave = false;
     }

     /*
         public void draw(GL gl, GLU glu, boolean isPlanar, float[] normal, StsColor color)
         {
             StsColor.setGLColor(gl, color);
             draw(gl, glu, isPlanar, normal);
         }
     */
     public void drawLines(GL gl)
     {
         gl.glBegin(GL.GL_LINE_LOOP);
         for (int n = 0; n < nPoints; n++)
             gl.glVertex3fv(pntNrmls[n], 0);
         gl.glEnd();
     }

     public String getLabel()
     {
         return new String("row: " + row + " col: " + col + " nPoints: " + nPoints);
     }

     public String toString()
     {
         float[] zRange = getZRange();
         int nPoints = pntNrmls != null ? pntNrmls.length : 0;
         return "row: " + row + " col: " + col + " layer: " + layer + " block: " + nBlock + " nPoints: " + nPoints +
             " type: " + typeStrings[type] + " zMin: " + zRange[0] + " zMax: " + zRange[1];
     }

     public float[] getZRange()
     {
         float zMin, zMax;
         if (pntNrmls != null)
         {
             zMin = pntNrmls[0][2];
             zMax = zMin;
             for (int n = 1; n < nPoints; n++)
             {
                 float z = pntNrmls[n][2];
                 if (z < zMin) zMin = z;
                 else if (z > zMax) zMax = z;
             }
         }
         else if (links != null && links.getSize() > 0)
         {
             StsEdgeLoopRadialGridLink link;
             link = (StsEdgeLoopRadialGridLink) links.getElement(0);
             zMin = link.getXYZ()[2];
             zMax = zMin;
             for (int n = 1; n < nPoints; n++)
             {
                 link = (StsEdgeLoopRadialGridLink) links.getElement(0);
                 float z = link.getXYZ()[2];
                 if (z < zMin) zMin = z;
                 else if (z > zMax) zMax = z;
             }
         }
         else
         {
             zMin = -99.0f;
             zMax = -99.0f;
         }
         return new float[]{zMin, zMax};
     }

     public int getID(){ return type; }

     public double[][] getPoints(){ return null; }

     public int getNPoints(){ return nPoints; }

     public int compareTo(StsPolygon other)
     {
         float[] zRange = getZRange();
         float[] otherZRange = other.getZRange();
         float minZ = zRange[0];
         float otherMinZ = otherZRange[0];
         if (minZ < otherMinZ) return -1;
         if (minZ > otherMinZ) return 1;
         float maxZ = zRange[1];
         float otherMaxZ = otherZRange[1];
         if (maxZ < otherMaxZ) return -1;
         if (maxZ > otherMaxZ) return 1;
         return 0;
     }

     private float getMinZ()
     {
         float minZ = StsParameters.largeFloat;
         for (int n = 0; n < nPoints; n++)
             if (pntNrmls[n][3] < minZ) minZ = pntNrmls[n][3];
         return minZ;
     }

     public float[][] getColSegment(int col, StsBlock block)
     {
         try
         {
             int nPoints = pntNrmls.length;
             int firstIndex = -1;
             int lastIndex = 0;
             for (int n = 0; n < nPoints; n++)
             {
                 if (block.isXonCol(pntNrmls[n][0], col))
                 {
                     if (firstIndex == -1)
                         firstIndex = n;
                     lastIndex = n;
                 }
             }
             if (firstIndex == -1) // no points on this col; return zero-length segment at nearest point
             {
                 float colX = block.getXCoor(col);
                 float distance = largeFloat;
                 for (int n = 0; n < nPoints; n++)
                 {
                     float xDistance = Math.abs(colX - pntNrmls[n][0]);
                     if (xDistance < distance)
                     {
                         firstIndex = n;
                         lastIndex = n;
                         distance = xDistance;
                     }
                 }
             }
             float[][] segment = new float[2][];
             segment[0] = pntNrmls[firstIndex];
             segment[1] = pntNrmls[lastIndex];
             return segment;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "getColSegment", e);
             return null;
         }
     }

     public float[][] getRowSegment(int row, StsBlock block)
     {
         int nPoints = pntNrmls.length;
         int firstIndex = -1;
         int lastIndex = -1;
         for (int n = 0; n <= nPoints; n++)
         {
             int index = n % nPoints;
             if (block.isYonRow(pntNrmls[index][1], row))
             {
                 if (firstIndex == -1)
                     firstIndex = index;
                 else
                     lastIndex = index;
             }
             else// if we only have one point in segment
             {
                 if (firstIndex != -1 && lastIndex == -1)
                     firstIndex = -1;
             }
         }
         if (firstIndex == -1 || lastIndex == -1) // no points on this col; return zero-length segment at nearest point
         {
             float rowY = block.getYCoor(row);
             float distance = largeFloat;
             for (int n = 0; n < nPoints; n++)
             {
                 float yDistance = Math.abs(rowY - pntNrmls[n][1]);
                 if (yDistance < distance)
                 {
                     firstIndex = n;
                     lastIndex = n;
                     distance = yDistance;
                 }
             }
         }
         float[][] segment = new float[2][];
         segment[0] = pntNrmls[firstIndex];
         segment[1] = pntNrmls[lastIndex];
         return segment;
     }

     static public float computePolygonFrustumVolumeOld(StsPolygon topPolygon, StsPolygon botPolygon)
     {
         float[][] topPoints = topPolygon.pntNrmls;
         float[][] botPoints = botPolygon.pntNrmls;
         int nSides = topPoints.length;
         if (nSides != botPoints.length)
         {
             StsException.systemError(StsPolygon.class, "computePolygonFrustVolume", "Unequal number of polygon sides.");
             return 0.0f;
         }
         double volume = 0.0;
         for (int n = 0; n < nSides; n++)
         {
             double tdx = topPoints[(n + 1) % nSides][0] - topPoints[n][0];
             double bdx = botPoints[(n + 1) % nSides][0] - botPoints[n][0];
             double tya = (topPoints[(n + 1) % nSides][1] + topPoints[n][1]) / 2;
             double bya = (botPoints[(n + 1) % nSides][1] + botPoints[n][1]) / 2;
             double tza = (topPoints[(n + 1) % nSides][2] + topPoints[n][2]) / 2;
             double bza = (botPoints[(n + 1) % nSides][2] + botPoints[n][2]) / 2;
             volume += (1.0 / 3.0) * (2 * (tdx * tya + bdx * bya) + tdx * bya + bdx * tya) * (bza - tza);
         }
         return (float) volume;
     }

     public double computePolygonArea()
     {
         float[][] points = pntNrmls;
         int nSides = points.length;
         double area = 0.0;
         float[] point = null;
         point = points[nSides - 1];
         for (int n = 0; n < nSides; n++)
         {
             float[] prevPoint = point;
             point = points[n];
             double tdx = point[0] - prevPoint[0];
             double tya = (point[1] + prevPoint[1]) / 2;
             area += tdx * tya;
         }
         if (area > 1.e10)
         {
             StsException.systemError(this, "computePolygonArea", "area is too large: " + area);
             return 0.0;
         }
         if (area < 0.0)
         {
             StsException.systemError(this, "computePolygonArea", "area is negative");
             return 0.0;
         }
         return area;
     }

     public double computeAverageZ()
     {
         int nSides = pntNrmls.length;
         double z = pntNrmls[0][2];
         for (int n = 1; n < nSides; n++)
             z += pntNrmls[n][2];
         return z / nSides;
     }

     public float[] computeCenter()
     {
         return StsMath.average(pntNrmls, 3);
     }

     float[] computePolygonCenterBad()
     {
         float[][] points = pntNrmls;
         float xc = 0, yc = 0, zc = 0;
         int nPoints = points.length;
         float[] point = points[nPoints - 1];
         float x = point[0];
         float y = point[1];
         float z = point[2];
         float zSum = 0.0f;
         float perimeter = 0.0f;
         for (int n = 0; n < nPoints; n++)
         {
             float xp = x;
             float yp = y;
             float zp = z;
             point = points[n];
             x = point[0];
             y = point[1];
             z = point[2];
             float dx = x - xp;
             float dy = y - yp;
             float sum = x * y + (dx * y + dy * x) / 2 + dx * dy / 3;
             xc += -dy * sum;
             yc += dx * sum;
             float length = (float) Math.sqrt(dx * dx + dy * dy);
             perimeter += length;
             zSum = length * (z + zp) / 2;
         }
         if (perimeter != 0.0f)
             zc = zSum / perimeter;
         else
             zc = points[0][2];
         return new float[]{xc, yc, zc};
     }

     static public float computePolygonVolume(StsPolygon polygon)
     {
         float[][] points = polygon.pntNrmls;
         int nSides = points.length;
         double volume = 0.0;
         float[] point = null;
         for (int n = nSides - 1; n >= 0; n--)
         {
             point = points[n];
             if (point != null) break;
         }
         if (point == null) return 0.0f;

         for (int n = 0; n < nSides; n++)
         {
             float[] prevPoint = point;
             while (points[n] == null)
             {
                 if (n >= nSides - 1) break;
                 n++;
             }
             if (points[n] == null) break;
             point = points[n];
             double tdx = point[0] - prevPoint[0];
             double tya = (point[1] + prevPoint[1]) / 2;
             double tza = (point[2] + prevPoint[2]) / 2;
             volume += tdx * tya * tza;
         }
         return (float) volume;
     }

     public int[] getIJKB(){ return new int[]{row, col, layer, nBlock}; }

     public double[] getNormal()
     {
         return new double[]{(double) pntNrmls[0][3], (double) pntNrmls[0][4], (double) pntNrmls[0][5]};
     }

     static public double computeIntersectionArea(StsPolygon polygonA, StsPolygon polygonB, double[] normal)
     {
         try
         {
             double[] Z = new double[]{0, 0, 1};
             double[] H = StsMath.crossNormalize(normal, Z);
             double[] V = StsMath.crossNormalize(normal, H);
             double[][] pA2 = StsPolygonIntersect.project2D(H, V, polygonA.pntNrmls);
             double[][] pB2 = StsPolygonIntersect.project2D(H, V, polygonB.pntNrmls);
             return StsPolygonIntersect.intersectionArea(pA2, pB2);
         }
         catch (Exception e)
         {
             StsException.systemError(StsPolygon.class, "computeIntersectionArea",
                 polygonA.toString() + " " + polygonB.toString() + StsMath.toString("normal", normal));
             return 0.0;
         }
     }

     static final int[][] shortDirections = new int[][]{{1, 2}, {0, 2}, {0, 1}};

     static public double computeIntersectionAreaNew(StsPolygon polygonA, StsPolygon polygonB, double[] normal)
     {
         try
         {
             int largestDirection = getLargestDirection(normal);
             double[][] pA2 = getPolygon2DFromDirection(shortDirections[largestDirection], polygonA.pntNrmls);
             double[][] pB2 = getPolygon2DFromDirection(shortDirections[largestDirection], polygonB.pntNrmls);

             double len = StsMath.length(normal); // length of normal

             return -StsPolygonIntersect.intersectionArea(pA2, pB2) * len / normal[largestDirection];
         }
         catch (Exception e)
         {
             StsException.systemError(StsPolygon.class, "computeIntersectionArea",
                 polygonA.toString() + " " + polygonB.toString() + StsMath.toString("normal", normal));
             return 0.0;
         }
     }

     static private int getLargestDirection(double[] normal)
     {
         double[] av = new double[3];
         for (int n = 0; n < 3; n++)
             av[n] = Math.abs(normal[n]);

         if (av[0] > av[1])
         {
             if (av[0] > av[2])                   // ignore x-coord

                 return 0;
         }
         else if (av[1] > av[2])                   // ignore y-coord

             return 1;

         return 2; // ignore z-coord
     }

     static private double[][] getPolygon2DFromDirection(int[] shortDirections, float[][] pntNrmls)
     {
         int nPoints = pntNrmls.length;
         double[][] points = new double[nPoints][2];
         int n0 = shortDirections[0];
         int n1 = shortDirections[1];
         for (int n = 0; n < nPoints; n++)
         {
             points[n][0] = pntNrmls[n][n0];
             points[n][1] = pntNrmls[n][n1];
         }
         return points;
     }

     public double computeArea2D()
     {
         float[][] points = pntNrmls;
         int nSides = points.length;
         double area = 0.0;
         float[] point = points[nSides - 2];
         float[] nextPoint = points[nSides - 1];
         for (int n = 0; n < nSides; n++)
         {
             float[] prevPoint = point;
             point = nextPoint;
             nextPoint = points[n];
             double x = point[0];
             double dy = (nextPoint[1] - prevPoint[1]);
             area -= x * dy;
         }
         if (area > 1.e10)
         {
             StsException.systemError(this, "computePolygonArea", "area is too large: " + area);
             return 0.0;
         }
         if (area < 0.0)
         {
             StsException.systemError(this, "computePolygonArea", "area is negative");
             return 0.0;
         }
         return area / 2;
     }

     public double computeArea2D(int n1, int n2)
     {
         float[][] points = pntNrmls;
         int nSides = points.length;
         double area = 0.0;
         float[] point = points[nSides - 2];
         float[] nextPoint = points[nSides - 1];
         for (int n = 0; n < nSides; n++)
         {
             float[] prevPoint = point;
             point = nextPoint;
             nextPoint = points[n];
             double x = point[n1];
             double dy = (nextPoint[n2] - prevPoint[n2]);
             area -= x * dy;
         }
         if (area > 1.e10)
         {
             StsException.systemError(this, "computePolygonArea", "area is too large: " + area);
             return 0.0;
         }
         if (area < 0.0)
         {
             StsException.systemError(this, "computePolygonArea", "area is negative");
             return 0.0;
         }
         return area / 2;
     }

     // return the signed area of a 3D planar polygon (given normal vector)

     public double computeArea3D(double[] normal)         // and plane normal

     {
         // select largest normal coordinate to ignore for projection

         double ax = Math.abs(normal[0]);
         double ay = Math.abs(normal[1]);
         double az = Math.abs(normal[2]);

         double len = StsMath.length(normal); // length of normal


         if (ax > ay)
         {
             if (ax > az)                   // ignore x-coord

                 return computeArea2D(1, 2) * (len / normal[0]);
         }
         else if (ay > az)                   // ignore y-coord

             return computeArea2D(0, 2) * (len / normal[1]);

         return computeArea2D(0, 1) * (len / normal[2]); // ignore z-coord
     }

     public StsColor setPropertyColor(StsPropertyVolume propertyVolume)
     {
         if (propertyVolume == null) return null;
         stsColor = propertyVolume.getColor(row, col, layer);
         return stsColor;
     }

     public void setGLColor(GL gl)
     {
         if (stsColor != null) stsColor.setGLColor(gl);
     }

     static public Iterator<StsPolygon> emptyIterator()
     {
         ArrayList<StsPolygon> emptyList = new ArrayList<StsPolygon>();
         return emptyList.iterator();
     }

     public ArrayList<StsPolygon> computeIntersectionPolygons(StsPolygon polygonB, double[] normal)
           {
               ArrayList<PolygonIntersectionPoint> intersectionPoints = computeIntersectionPoints(polygonB, normal);
               int nIntersectionPoints = intersectionPoints.size();
               if (nIntersectionPoints % 2 != 0)
                  StsException.systemError(this, "computeIntersectionPolygons", "number of intersections is odd.");

               PolygonIntersectionLoop loopA = new PolygonIntersectionLoop(intersectionPoints, true);
               PolygonIntersectionLoop loopB = new PolygonIntersectionLoop(intersectionPoints, false);

               ArrayList<StsPolygon> intersectionPolygons = new ArrayList<StsPolygon>();
               intersectionPoints = loopA.points;
               ArrayList<float[]> polygonPoints = new ArrayList<float[]>();
               for (int n = 0; n < nIntersectionPoints; n++)
               {
                   PolygonIntersectionPoint point = intersectionPoints.get(n);
                   while(!point.intersected)
                   {
                       StsPolygon polygon = new StsPolygon();
                        point.intersected = true;
                       polygonPoints.add(point.xyz);
                        // PolygonIntersectionPoint nextPoint = point.getNextPoint(loopA, loopB);

                   }
                   /*
                   PolygonIntersectionPoint pointA1 = pointA0.nextA;
                   if (pointA1.intersected) continue;
                   StsPolygon intersectionPolygon = constructIntersectionPolygon(polygonB, pointA0, pointA1);
                   if (intersectionPolygon != null)
                       intersectionPolygons.add(intersectionPolygon);
                       */
               }
               return intersectionPolygons;
           }

           private StsPolygon constructIntersectionPolygon(StsPolygon polygonB, PolygonIntersectionPoint point0, PolygonIntersectionPoint point1)
           {
               if (point0.inside == point1.inside)
               {
                   StsException.systemError(this, "constructIntersectionPolygon", "points are not opposing.");
                   return null;
               }
               if (point0.inside)
               {
                   int belowB1 = point1.getBelowB();
                   int aboveB0 = point0.getAboveB();
                   int nBPoints = belowB1 - aboveB0 + 1;
                   int nTotalAPoints = pntNrmls.length;
                   int belowA0 = nTotalAPoints + point0.getBelowA();
                   int aboveA1 = point1.getAboveA();
                   int nAPoints = belowA0 - aboveA1 + 1;
                   int nPoints = nAPoints + nBPoints + 2;
                   float[][] pntNrmls = new float[nPoints][];
                   int i = 0;
                   pntNrmls[i++] = point0.getPolygonPoint(this, polygonB);
                   for (int n = 0, p = aboveB0; n < nBPoints; n++, p++)
                       pntNrmls[i++] = this.pntNrmls[p];
                   pntNrmls[i++] = point1.getPolygonPoint(this, polygonB);
                   for (int n = 0, p = aboveA1; n < nAPoints; n++, p++)
                       pntNrmls[i++] = this.pntNrmls[p % nTotalAPoints];

                   StsPolygon polygon = new StsPolygon();
                   polygon.setPntNrmls(pntNrmls);
                   return polygon;
               }
               else // point1.inside
               {
                   int belowA1 = point1.getBelowA();
                   int aboveA0 = point0.getAboveA();
                   int nAPoints = belowA1 - aboveA0 + 1;
                   int belowB1 = point1.getBelowB();
                   int aboveB0 = point0.getAboveB();
                   int nBPoints = belowB1 - aboveB0 + 1;
                   int nPoints = nAPoints + nBPoints + 2;
                   float[][] pntNrmls = new float[nPoints][];
                   int i = 0;
                   pntNrmls[i++] = point0.getPolygonPoint(this, polygonB);
                   for (int n = 0, p = aboveA0; n < nAPoints; n++, p++)
                       pntNrmls[i++] = this.pntNrmls[p];
                   pntNrmls[i++] = point1.getPolygonPoint(this, polygonB);
                   for (int n = 0, p = belowB1; n < nBPoints; n++, p--)
                       pntNrmls[i++] = this.pntNrmls[p];
                   StsPolygon polygon = new StsPolygon();
                   polygon.setPntNrmls(pntNrmls);
                   return polygon;
               }
           }

           public ArrayList<PolygonIntersectionPoint> computeIntersectionPoints(StsPolygon polygonB, double[] normal)
           {
               ArrayList<PolygonIntersectionPoint> intersectionPoints = new ArrayList<PolygonIntersectionPoint>();
               float[][] pointsA = pntNrmls;
               float[][] pointsB = polygonB.pntNrmls;
               int largestDirection = getLargestDirection(normal);
               double[][] pA2 = getPolygon2DFromDirection(shortDirections[largestDirection], pointsA);
               double[][] pB2 = getPolygon2DFromDirection(shortDirections[largestDirection], pointsB);
               int nPointsA = pA2.length;
               int nPointsB = pB2.length;
               double[] factors = new double[2];
               double[] pointA1 = pA2[0];
               for (int a = 0; a < nPointsA; a++)
               {
                   double[] pointA0 = pointA1;
                   pointA1 = pA2[(a + 1) % nPointsA];

                   double[] pointB1 = pB2[0];
                   for (int b = 0; b < nPointsB; b++)
                   {
                       double[] pointB0 = pointB1;
                       pointB1 = pB2[(b + 1) % nPointsB];
                       PolygonIntersectionPoint intersectionPoint = checkAddIntersectionPoint(a, pointA0, pointA1, b, pointB0, pointB1, factors);
                       if (intersectionOk(intersectionPoint))
                           intersectionPoints.add(intersectionPoint);
                   }
               }
               return intersectionPoints;
           }

           private boolean intersectionOk(PolygonIntersectionPoint intersectionPoint)
           {
               if (intersectionPoint == null) return false;
               return intersectionPoint.isOk();
           }

           private PolygonIntersectionPoint checkAddIntersectionPoint(int nA, double[] pointA0, double[] pointA1, int nB, double[] pointB0, double[] pointB1, double[] factors)
           {
               double[] v0 = new double[]{pointA1[0] - pointA0[0], pointA1[1] - pointA0[1]};
               double[] v1 = new double[]{pointB1[0] - pointB0[0], pointB1[1] - pointB0[1]};
               double[] point = StsMath.linePVIntersectXY(pointA0, v0, pointB0, v1, factors);
               if (point == null) return null;
               boolean inside = StsMath.cross2(v1, v0) > 0.0;
               return new PolygonIntersectionPoint(nA, factors[0], nB, factors[1], inside);
           }

           class PolygonIntersectionPoint
           {
               float[] xyz;
               int nA, nB;
               double fA, fB;
               boolean inside;
               StsPolygon polygonA, polygonB;
               PolygonIntersectionPoint nextA, nextB;
               boolean intersected = false;

               PolygonIntersectionPoint(int nA, double fA, int nB, double fB, boolean inside)
               {
                   this.nA = nA;
                   this.fA = fA;
                   this.nB = nB;
                   this.fB = fB;
                   this.inside = inside;
               }

               boolean isOk()
               {
                   if (StsMath.betweenExclusive(fA, 0.0, 1.0) && StsMath.betweenExclusive(fB, 0.0, 1.0))
                       return true;
                   if (inside)
                       return fA >= 0.0 && fA < 1.0 && fB >= 0.0 && fB < 1.0;
                   else //outside
                       return fA > 0.0 && fA <= 1.0 && fB > 0.0 && fB <= 1.0;
               }

               int getBelowA()
               {
                   return StsMath.below(nA + fA);
               }

               int getBelowB()
               {
                   return StsMath.below(nB + fB);
               }

               int getAboveA()
               {
                   return StsMath.above(nA + fA);
               }

               int getAboveB()
               {
                   return StsMath.above(nB + fB);
               }

               float[] getPolygonPoint(StsPolygon polygonA, StsPolygon polygonB)
               {
                   if (fA == 0.0)
                       return polygonA.pntNrmls[nA];
                   if (fB == 0.0)
                       return polygonB.pntNrmls[nB];
                   if (fA == 1.0)
                   {
                       int n = (nA + 1) % (polygonA.pntNrmls.length);
                       return polygonA.pntNrmls[n];
                   }
                   if (fB == 1.0)
                   {
                       int n = (nB + 1) % (polygonB.pntNrmls.length);
                       return polygonB.pntNrmls[n];
                   }
                   else
                   {
                       float[][] pntNrmls = polygonA.pntNrmls;
                       int n0 = nA;
                       int n1 = (nA + 1) % (pntNrmls.length);
                       return StsMath.interpolate(pntNrmls[n0], pntNrmls[n1], (float) fA);
                   }
               }
           }

           class PolygonIntersectionLoop
           {
               ArrayList<PolygonIntersectionPoint> points;
               boolean isMainLoop;

               PolygonIntersectionLoop(ArrayList<PolygonIntersectionPoint> points, boolean isMainLoop)
               {
                   this.points = new ArrayList<PolygonIntersectionPoint>(points);
                   this.isMainLoop = isMainLoop;
                   sortPoints();
               }

               void sortPoints()
               {
                   int nPoints = points.size();
                   Collections.sort(points, getComparator());
                   PolygonIntersectionPoint nextPoint = points.get(0);
                   for (int n = 0; n < nPoints; n++)
                   {
                       PolygonIntersectionPoint prevPoint = nextPoint;
                       nextPoint = points.get((n + 1) % nPoints);
                       prevPoint.nextB = nextPoint;
                   }
               }

               private Comparator<PolygonIntersectionPoint> getComparator()
               {
                   if (isMainLoop)
                       return getComparatorA();
                   else
                       return getComparatorB();
               }

               private Comparator<PolygonIntersectionPoint> getComparatorA()
               {
                   return new Comparator<PolygonIntersectionPoint>()
                   {
                       public int compare(PolygonIntersectionPoint p1, PolygonIntersectionPoint p2)
                       {
                           return Double.compare(p1.nA + p1.fA, p2.nA + p2.fA);
                       }
                   };
               }

               private Comparator<PolygonIntersectionPoint> getComparatorB()
               {
                   return new Comparator<PolygonIntersectionPoint>()
                   {
                       public int compare(PolygonIntersectionPoint p1, PolygonIntersectionPoint p2)
                       {
                           return Double.compare(p1.nB + p1.fB, p2.nB + p2.fB);
                       }
                   };
               }
           }


           static public void main(String[] args)
           {
               float[][] pointsA, pointsB;
               StsPolygon polygonB;
               ArrayList<StsPolygon> polygons;

               // double[] normal = new double[] { 0.0, 0.0, 1.0 };
               double[] normal = new double[]{0.0, 0.0, 1.0};
               // points = new float[][] { { 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 1.0f, 0.0f, 0.0f } };
               pointsA = new float[][]{{0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {1.0f, 0.0f, 0.0f}};
               StsPolygon polygonA = initializeTestPolygon(pointsA, normal, "PolygonA");

               pointsB = new float[][]{{0.0f, 0.0f, 0.0f}, {0.0f, 0.5f, 0.0f}, {0.5f, 0.0f, 0.0f}};
               polygonB = initializeTestPolygon(pointsB, normal, "PolygonB");
               polygons = polygonA.computeIntersectionPolygons(polygonB, normal);
               processTestPolygonResults(0.125, normal, polygons);

               pointsB = new float[][]{{0.0f, 0.5f, 0.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {0.5f, 0.0f, 0.0f}};
               polygonB = initializeTestPolygon(pointsB, normal, "PolygonB");
               polygons = polygonA.computeIntersectionPolygons(polygonB, normal);
               processTestPolygonResults(0.875, normal, polygons);
           }

           static void processTestPolygonResults(double correctArea, double[] normal, ArrayList<StsPolygon> polygons)
           {
               for (StsPolygon polygon : polygons)
               {
                   System.out.println("Intersect polygon: " + polygon.toString());
                   double area = polygon.computeArea3D(normal);
                   System.out.println("area intersect(should be " + correctArea + "): " + area);
               }
           }

           static StsPolygon initializeTestPolygon(float[][] points, double[] normal, String name)
           {
               StsPolygon polygon = new StsPolygon();
               polygon.setPntNrmls(points);
               System.out.println("PolygonA: " + polygon.toString());
               double area = polygon.computeArea3D(normal);
               System.out.println(name + " " + +area);
               return polygon;
           }
       }
