//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
 import com.Sts.Interfaces.*;
 import com.Sts.MVC.View3d.*;
 import com.Sts.Utilities.*;

 import java.util.*;

public class StsEdgeLoopRadialGridLink
 {
     /** @param point point this link is associated with */
     protected StsGridSectionPoint point;
     /** @param edgeLoop loop this point belongs to */
     protected StsEdgeLoop edgeLoop;
     /** @param rowCol contains grid row & col info */
     protected StsRowCol rowCol;
     /** pointLinks in seq: ROW+, COL+, ROW-, COL-, (PREV, NEXT) */
     protected StsEdgeLoopRadialGridLink[] pointLinks;
     /** CCW sequence from PREV thru connected ROW/COL to NEXT */
     public byte[] linkedDirections;
     /** connectTypes of pointLinks */
     protected byte[] connectTypes;
     /** general purpose flag */
     protected byte flag;
     /** inside direction type; used by toString() */
     public byte insideDirection = NO_LINK;

     // State flags for connect
     static public final int MINUS = StsParameters.MINUS;
     static public final int PLUS = StsParameters.PLUS;

     // State flags for rowOrCol
     static final int ROW = StsParameters.ROW;
     static final int COL = StsParameters.COL;
     static final int ROWCOL = StsParameters.ROWCOL;
     static final int NONE = StsParameters.NONE;

     // flags for link directions
     static public final byte NO_LINK = StsParameters.NO_LINK;
     static public final byte MULTI_INSIDE_LINKS = StsParameters.MULTI_INSIDE_LINKS;

     static public final byte PREV_POINT = StsParameters.PREV_POINT;
     static public final byte NEXT_POINT = StsParameters.NEXT_POINT;

     static public final byte ROW_PLUS = StsParameters.ROW_PLUS;
     static public final byte COL_PLUS = StsParameters.COL_PLUS;
     static public final byte ROW_MINUS = StsParameters.ROW_MINUS;
     static public final byte COL_MINUS = StsParameters.COL_MINUS;

     static public final byte GRID_ROW_PLUS = StsParameters.GRID_ROW_PLUS;
     static public final byte GRID_COL_PLUS = StsParameters.GRID_COL_PLUS;
     static public final byte GRID_ROW_MINUS = StsParameters.GRID_ROW_MINUS;
     static public final byte GRID_COL_MINUS = StsParameters.GRID_COL_MINUS;

     static public final byte NOT_CONNECTED = 0;     // no connection in this direction
     static public final byte EDGE_CONNECTED = 1;
     static public final byte INSIDE_CONNECTED = 2;
     static public final byte OUTSIDE_CONNECTED = 3;
     static public final byte NEED_CONNECT = -1;
     static public final byte DONT_CONNECT = -2;


     static final float nullValue = StsParameters.nullValue;
     static final float roundOff = StsParameters.roundOff;
     static final float largeFloat = StsParameters.largeFloat;

     static public final byte ROW_USED = -126;
     static public final byte COL_USED = -127;
     static public final byte ROWCOL_USED = -128;

     public StsEdgeLoopRadialGridLink()
     {
     }

     /**
      * constructor for links on edges.  A total of six links possible: four directional links and two edge links.
      * For interior links use constructInteriorLink
      */
     public StsEdgeLoopRadialGridLink(StsEdgeLoop edgeLoop, StsGridSectionPoint gridPoint)
     {
         this(edgeLoop, gridPoint, 6, DONT_CONNECT);
     }

     public StsEdgeLoopRadialGridLink(StsEdgeLoop edgeLoop, StsGridSectionPoint gridPoint, int nLinks, byte connectType)
     {
         this.point = gridPoint;
         this.edgeLoop = edgeLoop;
         pointLinks = new StsEdgeLoopRadialGridLink[nLinks];
         initializeConnectTypes(nLinks, connectType);
         rowCol = gridPoint.getRowCol(edgeLoop.getGrid());
     }

     public StsEdgeLoopRadialGridLink(StsGridSectionPoint gridPoint, StsRowCol rowCol)
     {
         this(gridPoint, 4, rowCol);
     }

     public StsEdgeLoopRadialGridLink(StsGridSectionPoint gridPoint, int nLinks, StsRowCol rowCol)
     {
         this.point = gridPoint;
         pointLinks = new StsEdgeLoopRadialGridLink[nLinks];
         initializeConnectTypes(nLinks, DONT_CONNECT);
         this.rowCol = rowCol;
     }

     private StsEdgeLoopRadialGridLink(StsEdgeLoopRadialGridLink link)
     {
         this.point = link.point;
         this.edgeLoop = link.edgeLoop;
         this.rowCol = link.rowCol;
         this.flag = link.flag;
         this.insideDirection = link.insideDirection;

         StsEdgeLoopRadialGridLink[] pointLinks = link.pointLinks;
         if(pointLinks != null)
         {
             int nPointLinks = pointLinks.length;
             this.pointLinks = new StsEdgeLoopRadialGridLink[nPointLinks];
             System.arraycopy(pointLinks, 0, this.pointLinks, 0, nPointLinks);
         }
         byte[] linkedDirections = link.linkedDirections;
         if(linkedDirections != null)
         {
             int nLinkedDirections =linkedDirections.length;
             this.linkedDirections = new byte[nLinkedDirections];
             System.arraycopy(linkedDirections, 0, this.linkedDirections, 0, nLinkedDirections);
         }
         byte[] connectTypes = link.connectTypes;
         if(connectTypes != null)
         {
             int nConnectTypes = connectTypes.length;
             this.connectTypes = new byte[nConnectTypes];
             System.arraycopy(connectTypes, 0, this.connectTypes, 0, nConnectTypes);
         }

     }
     private void initializeConnectTypes(int nLinks, byte connectType)
     {
         connectTypes = new byte[nLinks];
         for (int n = 0; n < nLinks; n++) connectTypes[n] = connectType;
     }

     static public StsEdgeLoopRadialGridLink constructInteriorLink(StsEdgeLoop edgeLoop, StsGridSectionPoint gridPoint)
     {
         return new StsEdgeLoopRadialGridLink(edgeLoop, gridPoint, 4, NEED_CONNECT);
     }

     static public StsEdgeLoopRadialGridLink constructInteriorGridLink(StsEdgeLoop edgeLoop, StsGridSectionPoint gridPoint)
     {
         return new StsEdgeLoopRadialGridLink(edgeLoop, gridPoint, 10, NEED_CONNECT);
     }

     public void resetPoint(StsEdgeLoop edgeLoop, StsGridSectionPoint gridPoint)
     {
         this.point = gridPoint;
         this.edgeLoop = edgeLoop;
         rowCol = gridPoint.getRowCol(edgeLoop.getGrid());
     }

     public StsEdgeLoopRadialGridLink copy() { return new StsEdgeLoopRadialGridLink(this); }

     public StsGridSectionPoint getPoint(){ return point; }

     public StsEdgeLoopRadialGridLink[] getPointLinks(){ return pointLinks; }

     public StsRowCol getRowCol(){ return rowCol; }

     public int getRowOrCol(){ return rowCol.getRowOrCol(); }

     public float getRowF(){ return rowCol.getRowF(); }

     public float getColF(){ return rowCol.getColF(); }

     public float getGridRowF(){ return point.getGridRowF(); }

     public float getGridColF(){ return point.getGridColF(); }

     // public int getID() { return edgeLoop.getIndex(); }
     public StsEdgeLoop getEdgeLoop()
     { return edgeLoop; }

     public byte getFlag(){ return flag; }

     public void setFlag(byte flag)
     {
         this.flag = flag;
     }

     public void setRowColUsedFlag(int rowOrCol)
     {
         if (!rowColUseOK(rowOrCol))
         {
             StsException.systemError(this, "setRowColUsedFlag", "Can't set flag to " + flagString(flag) +
                 " already set to " + flagString(this.flag));
             return;
         }
         if (rowOrCol == ROW)
         {
             if (this.flag == COL_USED)
                 this.flag = ROWCOL_USED;
             else
                 this.flag = ROW_USED;
         }
         else // rowOrCol == COL
         {
             if (this.flag == ROW_USED)
                 this.flag = ROWCOL_USED;
             else
                 this.flag = COL_USED;
         }
     }

     static public final String flagString(byte flag)
     {
         switch (flag)
         {
             case ROW_USED:
                 return "ROW_USED";
             case COL_USED:
                 return "COL_USED";
             case ROWCOL_USED:
                 return "ROWCOL_USED";
             default:
                 return "UNKNOWN";
         }
     }

     public boolean rowColUsedFlagOK(byte newFlag)
     {
         if (this.flag == ROWCOL_USED)
             return false;
         else
             return this.flag != newFlag;
     }

     public boolean rowColUseOK(int rowOrCol)
     {
         if (this.flag == rowOrCol) return true;

         if (this.flag == ROWCOL_USED)
             return false;
         if (rowOrCol == ROW)
             return this.flag != ROW_USED;
         else  // rowOrCol == COL
             return this.flag != COL_USED;
     }

     public boolean isRowOrCol()
     {
         return rowCol.getRowOrCol() != NONE;
     }

     public boolean isRowAndCol()
     {
         return rowCol.getRowOrCol() == ROWCOL;
     }

     public boolean isRow()
     {
         int rowOrCol = rowCol.getRowOrCol();
         return rowOrCol == ROW || rowOrCol == ROWCOL;
     }

     public boolean isCol()
     {
         int rowOrCol = rowCol.getRowOrCol();
         return rowOrCol == COL || rowOrCol == ROWCOL;
     }

     public int getRow()
     {
         float rowF = rowCol.getRowF();
         return StsMath.roundOffInteger(rowF);
     }

     public int getRowOrColNumber(int rowOrCol)
     {
         if (rowOrCol == ROW)
             return getRow();
         else if (rowOrCol == COL)
             return getCol();
         StsException.systemError(this, "getRowOrColNumber", "Must be called with 'ROW' or 'COL', but called with " + StsParameters.rowCol(rowOrCol));
         return -1;
     }

     public int getCol()
     {
         float colF = rowCol.getColF();
         return StsMath.roundOffInteger(colF);
     }

     static public int getLinkIndexRowOrCol(int linkIndex)
     {
         if (linkIndex == 0 || linkIndex == 2)
             return ROW;
         else if (linkIndex == 1 || linkIndex == 3)
             return COL;
         else
             return NONE;
     }

     static public int getLinkIndexDirection(int linkIndex)
     {
         if (linkIndex == 0 || linkIndex == 1)
             return PLUS;
         else if (linkIndex == 2 || linkIndex == 3)
             return MINUS;
         else
             return NONE;
     }

     static public int getOppositeLinkIndex(int linkIndex) throws StsException
     {
         switch (linkIndex)
         {
             case 0:
                 return 2;
             case 1:
                 return 3;
             case 2:
                 return 0;
             case 3:
                 return 1;
             case 6:
                 return 8;
             case 7:
                 return 9;
             case 8:
                 return 6;
             case 9:
                 return 7;
             default:
                 throw new StsException(StsException.WARNING, "StsEdgeLoopRadialGridLink.getOppositeLinkIndex() failed.",
                     " Incorrect input linkIndex: " + linkIndex);
         }
     }

     public float[] getXYZ()
     {
         return point.getPoint().getXYZorT();
     }

     public float[] getNormal(StsSurfaceGridable grid)
     {
         return grid.getNormal(getRowF(), getColF());
     }

     public float[] getVertexAndNormal(StsSurfaceGridable grid)
     {
         try
         {
             float[] vertex = point.getPoint().getXYZorT();
             float[] normal = grid.getNormal(getRowF(), getColF());
             if (normal != null)
                 return new float[]{vertex[0], vertex[1], vertex[2], normal[0], normal[1], normal[2]};
             else
                 return new float[]{vertex[0], vertex[1], vertex[2], 0.0f, 0.0f, 1.0f};
         }
         catch (Exception e)
         {
             StsException.systemError("StsEdgeLoopRadialGridLink.getVertexAndNormal() failed." + toString());
             return null;
         }
     }

     public boolean isRowOrCol(int matchRowOrCol)
     {
         int rowOrCol = rowCol.getRowOrCol();
         if (matchRowOrCol == rowOrCol) return true;

         if (matchRowOrCol == ROW)
             return rowOrCol == ROWCOL;
         else if (matchRowOrCol == COL)
             return rowOrCol == ROWCOL;
         else
             return false;
     }

     public boolean isSameRowAndCol(StsEdgeLoopRadialGridLink otherLink)
     {

         if (otherLink == null) return false;
         StsRowCol otherRowCol = otherLink.getRowCol();

         if (rowCol.getRowF() != otherRowCol.getRowF()) return false;
         return rowCol.getColF() == otherRowCol.getColF();
     }

     public boolean isOnSameRowOrCol(StsEdgeLoopRadialGridLink otherLink)
     {

         float dRow = getRowF() - otherLink.getRowF();
         if (dRow == 0.0f && rowCol.isRowOrCol(ROW)) return true;

         float dCol = getColF() - otherLink.getColF();
         if (dCol == 0.0f && rowCol.isRowOrCol(COL)) return true;

         return false;
     }

     public boolean isSameRowAndCol(int row, int col)
     {
         if (rowCol.getRowF() != (float) row) return false;
         return rowCol.getColF() == (float) col;
     }

     public boolean isRowOrColAligned(StsEdgeLoopRadialGridLink otherLink)
     {
         if (otherLink == null) return false;
         StsRowCol otherRowCol = otherLink.getRowCol();

         return rowCol.getRowF() == otherRowCol.getRowF() ||
             rowCol.getColF() == otherRowCol.getColF();
     }

     public boolean isSameRow(StsSurfaceGridable grid, StsEdgeLoopRadialGridLink otherLink)
     {
         if (otherLink == null) return false;
         float rowF = point.getRowF(grid);
         float otherRowF = otherLink.getPoint().getRowF(grid);
         return rowF == otherRowF;
     }

     public boolean isSameCol(StsSurfaceGridable grid, StsEdgeLoopRadialGridLink otherLink)
     {
         if (otherLink == null) return false;
         float colF = point.getColF(grid);
         float otherColF = otherLink.getPoint().getColF(grid);
         return colF == otherColF;
     }

     public float getRowColF(int rowOrCol)
     {
         if (rowOrCol == ROW) return rowCol.getRowF();
         else if (rowOrCol == COL) return rowCol.getColF();
         else return nullValue;
     }

     public int getRowCol(int rowOrCol)
     {
         if (rowOrCol == ROW) return (int) rowCol.getRowF();
         else if (rowOrCol == COL) return (int) rowCol.getColF();
         else return -1;
     }

     public int getGridRowCol(int rowOrCol)
     {
         if (rowOrCol == ROW) return (int) getGridRowF();
         else if (rowOrCol == COL) return (int) getGridColF();
         else return -1;
     }

     public float getCrossingRowColF()
     {
         return rowCol.getCrossingRowColF();
     }

     public float getCrossingRowColF(int rowOrCol)
     {
         return rowCol.getCrossingRowColF(rowOrCol);
     }

     public boolean matchesRowCol(int[] rowColIndexes)
     {
         return rowCol.matchesRowCol(rowColIndexes);
     }

     public boolean isSameEdgeLoop(StsEdgeLoopRadialGridLink otherLink)
     {
         return edgeLoop == otherLink.getEdgeLoop();
     }

     /*
      public void setRowOrCol(int setRowOrCol, float rowCol)
         {
             if(rowOrCol == NONE)
                 rowOrCol = setRowOrCol;
             else if(rowOrCol != setRowOrCol)
                 rowOrCol = ROWCOL;

             if(setRowOrCol == ROW)
                 rowF = rowCol;
             else if(setRowOrCol == COL)
                 colF = rowCol;
         }

         public void computeRowOrCol()
         {
             boolean isRow = StsMath.isIntegral(rowF);
             boolean isCol = StsMath.isIntegral(colF);

             if(isRow) rowF = StsMath.roundOffInteger(rowF);
             if(isCol) colF = StsMath.roundOffInteger(colF);

             if(isRow)
             {
                 if(isCol) rowOrCol = ROWCOL;
                 else      rowOrCol = ROW;
             }
             else if(isCol)
                 rowOrCol = COL;
             else
                 rowOrCol = NONE;
         }

         public void setRowColF(float rowF, float colF)
         {
             this.rowF = rowF;
             this.colF = colF;
             computeRowOrCol();
         }

         public void setRowColF(int rowOrCol, float rowF, float colF)
         {
             this.rowOrCol = rowOrCol;
             this.rowF = rowF;
             this.colF = colF;
         }
     */
     public boolean isConnected(int linkIndex)
     {
         if (linkIndex == NO_LINK) return false;
         return connectTypes[linkIndex] > 0;
     }

     public boolean isConnected(int rowOrCol, int direction)
     {
         return isConnected(getLinkDirection(rowOrCol, direction));
     }

     public boolean needsConnect(int rowOrCol)
     {
         if (rowOrCol == ROW &&
             (connectTypes[ROW_MINUS] == NEED_CONNECT || connectTypes[ROW_PLUS] == NEED_CONNECT))
             return true;

         if (rowOrCol == COL &&
             (connectTypes[COL_MINUS] == NEED_CONNECT || connectTypes[COL_PLUS] == NEED_CONNECT))
             return true;
         return false;
     }

     // if there is not exactly one inside direction, we have none or multiple: return NONE
     public byte getInsideDirection()
     {
        return insideDirection;
     }

     // Check all the inside directions to see which correspond to this rowOrCol.
     // If there is exactly 1, return the corresponding direction, otherwise return NONE;
     public int getInsideDirection(int rowOrCol)
     {
         int outsideDirection = getOutsideDirection(rowOrCol);
         if (outsideDirection == NONE) return NONE;
         else return -outsideDirection;
     }

     /*
      public int getInsideDirection(int rowOrCol)
         {
             int[] insideLinkIndexes = getInsideLinkIndexes();
             boolean isPlus = false;
             boolean isMinus = false;

             for(int n = 0; n < insideLinkIndexes.length; n++)
             {
                 int direction = getDirection(insideLinkIndexes[n], rowOrCol);
                 if(direction == PLUS) isPlus = true;
                 else if(direction == MINUS) isMinus = true;
             }

             if(isPlus && isMinus || !isPlus && !isMinus) return NONE;
             else if(isPlus) return PLUS;
             else return MINUS;
         }
     */
     private int getDirection(int linkIndex, int rowOrCol)
     {
         if (rowOrCol == ROW)
         {
             if (linkIndex == 0) return PLUS;
             else if (linkIndex == 2) return MINUS;
             else return NONE;
         }
         else if (rowOrCol == COL)
         {
             if (linkIndex == 1) return PLUS;
             else if (linkIndex == 3) return MINUS;
             else return NONE;
         }
         else if (rowOrCol == ROWCOL)
         {
             if (linkIndex == 0 || linkIndex == 1) return PLUS;
             else return MINUS;
         }
         else
             return NONE;
     }

     public int getOutsideDirection(int rowOrCol)
     {
         boolean minusNotConnected = getConnectType(rowOrCol, MINUS) == DONT_CONNECT;
         boolean plusNotConnected = getConnectType(rowOrCol, PLUS) == DONT_CONNECT;

         if (plusNotConnected && minusNotConnected || !plusNotConnected && !minusNotConnected)
             return NONE;
         else if (plusNotConnected)
             return PLUS;
         else
             return MINUS;
     }

     public boolean hasOutsideDirection(int rowOrCol, int direction)
     {
         int outsideDirection = getOutsideDirection(rowOrCol);
         return direction == outsideDirection;
     }
     /*
         public int getOutsideDirection(int rowOrCol)
            {
                int insideDirection = getInsideDirection(rowOrCol);
                if(insideDirection == NONE) return NONE;
                else return -insideDirection;
            }
        */

     public boolean hasInsideDirection(int rowOrCol, int direction)
     {
         int linkIndex = StsEdgeLoopRadialGridLink.getLinkDirection(rowOrCol, direction);
         return hasInsideDirection(linkIndex);
     }

     public boolean hasInsideDirection(int insideDirection)
     {
         if(linkedDirections == null || linkedDirections.length == 0)
            return false;
         for (int n = 0; n < linkedDirections.length; n++)
             if (linkedDirections[n] == insideDirection) return true;
         return false;
     }

     public StsEdgeLoopRadialGridLink getNextEdgeLink()
     {
         if (pointLinks == null || pointLinks.length < 6)
             return null;
         return pointLinks[NEXT_POINT];
     }

     public StsEdgeLoopRadialGridLink getPrevEdgeLink()
     {
         if (pointLinks == null || pointLinks.length < 6)
             return null;
         return pointLinks[PREV_POINT];
     }

     public boolean isLinked(int rowOrCol, int direction)
     {
         return isLinked(getLinkDirection(rowOrCol, direction));
     }

     public boolean isGridLinked(int rowOrCol, int direction)
     {
         return isLinked(getGridLinkIndex(rowOrCol, direction));
     }

     public boolean isLinked(int linkIndex)
     {
         if (linkIndex == NO_LINK) return false;
         return pointLinks[linkIndex] != null;
     }

     /*
      public void initializeLinks(StsEdgeLoopRadialGridLink prevLink, StsEdgeLoopRadialGridLink nextLink)
         {
             int n, i;

             if(pointLinks == null || pointLinks.length < 6)
             {
                 pointLinks = new StsEdgeLoopRadialGridLink[6];
                 connectTypes = new byte[6];
                 for(n = 0; n < 6; n++) connectTypes[n] = DONT_CONNECT;
             }

             if(prevLink != null) pointLinks[PREV_POINT] = prevLink;
             if(nextLink != null) pointLinks[NEXT_POINT] = nextLink;
         }
     */
     public void setPrevLink(StsEdgeLoopRadialGridLink prevLink)
     {
         pointLinks[PREV_POINT] = prevLink;
     }

     public void setNextLink(StsEdgeLoopRadialGridLink nextLink)
     {
         pointLinks[NEXT_POINT] = nextLink;
     }

     public void addTwoWayGridConnects()
     {
         int n, i;

         int rowOrCol = rowCol.getRowOrCol();
         if (rowOrCol == NONE) return;

         switch (rowOrCol)
         {
             case ROWCOL:
                 for (n = 0; n < 4; n++) connectTypes[n] = NEED_CONNECT;
                 break;
             case ROW:
                 connectTypes[0] = NEED_CONNECT;
                 connectTypes[2] = NEED_CONNECT;
                 break;
             case COL:
                 connectTypes[1] = NEED_CONNECT;
                 connectTypes[3] = NEED_CONNECT;
                 break;
         }
     }

     public void addGridConnects()
     {
         int n, i;
         int connectStartIndex, connectEndIndex;
         int linkStartIndex, linkEndIndex;

         try
         {
             if (rowCol.getRowOrCol() == NONE) return;

             StsEdgeLoopRadialGridLink prevLink = pointLinks[PREV_POINT];
             StsEdgeLoopRadialGridLink nextLink = pointLinks[NEXT_POINT];
             if (prevLink == null || nextLink == null) return;

             connectTypes[PREV_POINT] = EDGE_CONNECTED;
             connectTypes[NEXT_POINT] = EDGE_CONNECTED;

             float prevAngle = getQuadAngle(prevLink);
             float nextAngle = getQuadAngle(nextLink);

             //        edgeAngles = new float[] { prevAngle, nextAngle };

             connectStartIndex = StsMath.ceiling(prevAngle);
             connectEndIndex = StsMath.floor(nextAngle);
             linkStartIndex = StsMath.above(prevAngle);
             linkEndIndex = StsMath.below(nextAngle);

             if (nextAngle < prevAngle)
             {
                 connectEndIndex += 4;
                 linkEndIndex += 4;
             }

             for (i = connectStartIndex; i <= connectEndIndex; i++)
             {
                 int direction = i % 4;
                 if (directionMatchesRowOrCol(direction)) connectTypes[direction] = NEED_CONNECT;
             }

             byte[] linkedDirections = new byte[6];

             int nLinks = 0;
             linkedDirections[nLinks++] = PREV_POINT;
             for (i = linkStartIndex; i <= linkEndIndex; i++)
             {
                 //                int direction = (i+4)%4;
                 int direction = i % 4;
                 if (direction < 0 || direction >= connectTypes.length)
                 {
                     StsException.systemError("StsEdgeLoopRadialGridLink error: direction " + direction + " outside connectTypes range (0 to " + connectTypes.length + ")");
                     continue;
                 }
                 if (connectTypes[direction] == NEED_CONNECT) linkedDirections[nLinks++] = (byte) direction;
             }
             linkedDirections[nLinks++] = NEXT_POINT;

             setLinkedDirections(linkedDirections, nLinks);
         }
         catch (Exception e)
         {
             StsException.outputException("StsEdgeLoopRadialGridLink.addGridConnects() failed.",
                 e, StsException.WARNING);
         }
     }
     /*
         public void addGridConnects()
            {
                int n, i;
                int connectStartIndex, connectEndIndex;
                int linkStartIndex, linkEndIndex;

                try
                {
                    if(rowCol.getRowOrCol() == NONE) return;

                    StsEdgeLoopRadialGridLink prevLink = pointLinks[PREV_POINT];
                    StsEdgeLoopRadialGridLink nextLink = pointLinks[NEXT_POINT];
                    if(prevLink == null || nextLink == null) return;

                    connectTypes[PREV_POINT] = EDGE_CONNECTED;
                    connectTypes[NEXT_POINT] = EDGE_CONNECTED;

                    float prevAngle = getQuadAngle(prevLink);
                    float nextAngle = getQuadAngle(nextLink);

            //        edgeAngles = new float[] { prevAngle, nextAngle };

                    connectStartIndex = StsMath.ceiling(prevAngle);
                    connectEndIndex = StsMath.floor(nextAngle);
                    linkStartIndex = StsMath.above(prevAngle);
                    linkEndIndex = StsMath.below(nextAngle);

                    if(nextAngle < prevAngle)
                    {
                        connectEndIndex += 4;
                        linkEndIndex += 4;
                    }

                    for(i = connectStartIndex; i <= connectEndIndex; i++)
                    {
                        int direction = i%4;
                        if(directionMatchesRowOrCol(direction)) connectTypes[direction] = NEED_CONNECT;
                    }

                    setLinkedDirections(new byte[] { PREV_POINT, NEXT_POINT });

                    byte[] tempLinkedDirections = new byte[6];

                    int nLinks = 0;
                    tempLinkedDirections[nLinks++] = PREV_POINT;
                    for(i = linkStartIndex; i <= linkEndIndex; i++)
                    {
        //                int direction = (i+4)%4;
                        int direction = i%4;
                        if(direction < 0 || direction >= connectTypes.length)
                        {
                            StsException.systemError("StsEdgeLoopRadialGridLink error: direction " + direction + " outside connectTypes range (0 to " + connectTypes.length + ")");
                            continue;
                        }
                        if(connectTypes[direction] == NEED_CONNECT) tempLinkedDirections[nLinks++] = (byte)direction;
                    }
                    tempLinkedDirections[nLinks++] = NEXT_POINT;

                    setLinkedDirections(new byte[nLinks]);
                    System.arraycopy(tempLinkedDirections, 0, getLinkedDirections(), 0, nLinks);
                }
                catch(Exception e)
                {
                    StsException.outputException("StsEdgeLoopRadialGridLink.addGridConnects() failed.",
                        e, StsException.WARNING);
                }
            }
        */

     // return the angle in degrees to the otherLink divided by 90 which is the quad it is in:
     // 0 - 1 for NE, 1 - 2 for NW, 2 - 3 for SW, 3 - 4 for SE

     private float getQuadAngle(StsEdgeLoopRadialGridLink otherLink)
     {
         int sign;

         try
         {
             StsRowCol otherRowCol = otherLink.getRowCol();

             float dRow = otherRowCol.getRowF() - getRowF();
             float dCol = otherRowCol.getColF() - getColF();

             if (dRow != 0.0f && dCol != 0.0f)
                 return (float) StsMath.atan2(dCol, dRow) / 90;
             else if (dRow == 0.0f)
             {
                 if (dCol == 0.0f) return -1.0f;
                 else if (dCol > 0.0f) return 0.0f;
                 else return 2.0f;
             }
             else // dCol == 0.0f && dRow != 0.0f
             {
                 if (dRow > 0.0f) return 1.0f;
                 else return 3.0f;
             }
         }
         catch (Exception e)
         {
             StsException.outputException("StsEdgeLoopRadialGridLink.getQuadAngle() failed.",
                 e, StsException.WARNING);
             return 0.0f;
         }
     }

     public boolean directionMatchesRowOrCol(int linkIndex)
     {
         int rowOrCol = rowCol.getRowOrCol();

         if ((rowOrCol == ROW || rowOrCol == ROWCOL) &&
             (linkIndex == ROW_MINUS || linkIndex == ROW_PLUS)) return true;

         if ((rowOrCol == COL || rowOrCol == ROWCOL) &&
             (linkIndex == COL_MINUS || linkIndex == COL_PLUS)) return true;

         return false;
     }

     public byte[] getLinkedDirections()
     {
         if (!checkConstructLinkedDirections(false)) return new byte[0];
         if (linkedDirections.length > 0)
             return linkedDirections;
         else
             return new byte[0];
     }

     public int[] getInsideLinkIndexesNew()
     {
         if (!checkConstructLinkedDirections(true)) return null;
         //        if(linkedDirections == null && !constructLinkedDirections()) return new int[0];
         int nLinks = linkedDirections.length;
         if (nLinks <= 0) return null;
         if (linkedDirections[0] != 4 || linkedDirections[nLinks - 1] != 5) return null;
         int nInsideLinks = nLinks - 2;
         int[] insideLinks = new int[nInsideLinks];
         for (int n = 1; n < nLinks - 1; n++)
             insideLinks[n - 1] = n;
         return insideLinks;
     }

     // there might be two insideGridLinks if this is a cusp, but return only the first
     public StsEdgeLoopRadialGridLink getInsideGridLink(int rowOrCol)
     {
         if (pointLinks.length < 10) return null;
         if (rowOrCol == ROW || rowOrCol == ROWCOL)
         {
             if (pointLinks[6] != null) return pointLinks[6];
             if (pointLinks[8] != null) return pointLinks[8];
             return null;
         }
         else if (rowOrCol == COL || rowOrCol == ROWCOL)
         {
             if (pointLinks[7] != null) return pointLinks[7];
             if (pointLinks[9] != null) return pointLinks[9];
             return null;
         }
         else
         {
             StsException.systemError(this, "getInsideGridLink", "called with bad rowOrCol: " + rowOrCol);
             return null;
         }
     }

     public StsEdgeLoopRadialGridLink getInsideIndexedGridLink(int index)
     {
         if (pointLinks == null || pointLinks.length <= index) return null;
         StsEdgeLoopRadialGridLink link = pointLinks[index];
         if (link == getPrevEdgeLink()) return null;
         if (link == getNextEdgeLink()) return null;
         return link;
     }

     public void addEdgeGridRowAndColLinks()
     {
         int nLinks = pointLinks.length;
         if (nLinks > 6) return;
         StsGridRowCol gridRowCol = point.getGridRowCol();
         boolean isRow = gridRowCol.isRow();
         boolean isCol = gridRowCol.isCol();
         if (!isRow && !isCol) return;
         StsEdgeLoopRadialGridLink nextLink = getNextEdgeLink();
         StsGridRowCol nextGridRowCol = nextLink.point.getGridRowCol();
         if (isRow && nextGridRowCol.isRow())
         {
             float rowF = gridRowCol.rowF;
             float otherRowF = nextGridRowCol.rowF;
             if (rowF == otherRowF)
             {
                 float colF = gridRowCol.getColF();
                 float otherColF = nextGridRowCol.getColF();
                 if (otherColF > colF)
                     connectLinks(nextLink, GRID_ROW_PLUS);
                 else if (colF > otherColF)
                     connectLinks(nextLink, GRID_ROW_MINUS);
             }
         }
         if (isCol && nextGridRowCol.isCol())
         {
             float colF = gridRowCol.colF;
             float otherColF = nextGridRowCol.colF;
             if (colF == otherColF)
             {
                 float rowF = gridRowCol.getRowF();
                 float otherRowF = nextGridRowCol.getRowF();
                 if (otherRowF > rowF)
                     connectLinks(nextLink, GRID_COL_PLUS);
                 else if (rowF > otherRowF)
                     connectLinks(nextLink, GRID_COL_MINUS);
             }
         }
         else
         {

         }
     }

     public boolean connectLinks(StsEdgeLoopRadialGridLink otherLink, int rowColIndex)
     {
         try
         {
             if (!connectLink(otherLink, rowColIndex)) return false;
             rowColIndex = StsEdgeLoopRadialGridLink.getOppositeLinkIndex(rowColIndex);
             if (!otherLink.connectLink(this, rowColIndex)) return false;
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("{StsEdgeLoop.connectLinks() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     public boolean connectLink(StsEdgeLoopRadialGridLink otherLink, int rowColIndex)
     {
         try
         {
             StsEdgeLoopRadialGridLink currentLink = getLink(rowColIndex);
             if (currentLink != null)
             {
                 if (this == otherLink) return true;
                 rowColIndex = StsEdgeLoopRadialGridLink.getOppositeLinkIndex(rowColIndex);
                 currentLink = getLink(rowColIndex);
                 if (currentLink != null)
                 {
                     StsException.systemError(this, "connectLinks", toString() + " cannot connect to: " + otherLink.toString() +
                         ". Already connected to: " + currentLink.toString());
                     return false;
                 }
             }
             setLink(rowColIndex, otherLink);
             setConnectType(rowColIndex, StsEdgeLoopRadialGridLink.INSIDE_CONNECTED);
             constructLinkedDirections();
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("{StsEdgeLoop.connectLinks() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }


     /*
 public int getNInsideLinks()
         {
             if(linkedDirections == null) return 0;
             return linkedDirections.length - 2;
         }

         public int getFirstOutsideLinkIndex()
         {
             if(linkedDirections == null) return NO_LINK;
             int nLinkedDirections = linkedDirections.length;
             int nInsideLinks = nLinkedDirections-2;
             if(nInsideLinks <= 0) return NO_LINK;

             for(int n = 1; n < nLinkedDirections-1; n++)
             {
                 int insideLink = linkedDirections[n];
                 if(insideLink == ROW_PLUS || insideLink == COL_PLUS) return insideLink;
             }
             return NO_LINK;
         }

         public int getFirstInsideLinkIndex()
         {
             if(linkedDirections == null) return -1;
             int nLinkedDirections = linkedDirections.length;
             int nInsideLinks = nLinkedDirections-2;
             if(nInsideLinks <= 0) return -1;
             return linkedDirections[1];
         }
     */
     public boolean hasInsideLink(StsEdgeLoopRadialGridLink prevEdgeLink, StsEdgeLoopRadialGridLink nextEdgeLink)
     {
         StsEdgeLoopRadialGridLink insideLink = getNextCounterClockwiseLink(prevEdgeLink);
         return insideLink != null && insideLink != nextEdgeLink;
     }

     /*
 /** for this link, determine the number of CCW links between the prev and next edgeLinks.  These are the number of inside connections. */
     public int getNumberOfInsideLinks(StsEdgeLoopRadialGridLink prevEdgeLink, StsEdgeLoopRadialGridLink nextEdgeLink)
     {
         StsEdgeLoopRadialGridLink link = getNextCounterClockwiseLink(prevEdgeLink);
         if (link == null) return 0;

         int nInsideConnections = 0;
         while (link != null && link != nextEdgeLink && link != prevEdgeLink)
         {
             nInsideConnections++;
             link = getNextCounterClockwiseLink(link);
         }
         return nInsideConnections;
     }

     public int getNumberOfLinkConnectionsNew()
     {
         if (linkedDirections == null) return 0;
         int nLinkedDirections = linkedDirections.length;
         if (nLinkedDirections == 0) return 0;
         if (linkedDirections[0] != 4 && linkedDirections[nLinkedDirections - 1] != 5) return 0;
         int nInsideConnections = linkedDirections.length - 2;
         nInsideConnections = Math.min(2, nInsideConnections);
         return 2 + 2 * nInsideConnections;
     }
     /*
 public boolean XhasInsideLinks(StsEdgeLoopRadialGridLink prevEdgeLink, StsEdgeLoopRadialGridLink nextEdgeLink)
         {
             int n, nn;

             if(pointLinks == null) return false;
             if(linkedDirections != null) return hasInsideLinks();

             // There are only 4 point links on row/col plus/minus.
             // Find indexes of prevEdgeLink and nextEdgeLink and return true if there are
             // links in between.

             int prevIndex = -1, nextIndex = -1;
             for(n = 0; n < 4; n++)
             {
                 if(pointLinks[n] == prevEdgeLink) prevIndex = n;
                 if(pointLinks[n] == nextEdgeLink) nextIndex = n;
             }

             // between these two, find the number of non-null edges: if > 0 return true
             if(prevIndex == -1 || nextIndex == -1) return false;
             if(prevIndex == nextIndex) return false;

             for(n = 0, nn = prevIndex+1; n < 4; n++, nn++)
             {
                 int i = nn%4;
                 if(i == nextIndex) return false;
                 if(pointLinks[i] != null) return true;
             }
             return false;
         }
     */

     /*
 public boolean hasInsideLinks()
         {
             if(linkedDirections == null && !constructLinkedDirections()) return false;
             int nLinkedDirections = linkedDirections.length;
             return nLinkedDirections > 2;
         }
     */

     public int getLinkIndex()
     {
         int rowOrCol = rowCol.getRowOrCol();
         int direction = getInsideDirection(rowOrCol);
         return getLinkDirection(rowOrCol, direction);
     }

     static public byte getLinkDirection(int rowOrCol, int direction)
     {
         if (direction == NONE) return NO_LINK;

         if (rowOrCol == NONE)
         {
             if (direction == PLUS) return NEXT_POINT;
             else return PREV_POINT;
         }
         else if (rowOrCol == ROW)
         {
             if (direction == PLUS) return ROW_PLUS;
             else return ROW_MINUS;
         }
         else if (rowOrCol == COL)
         {
             if (direction == PLUS) return COL_PLUS;
             else return COL_MINUS;
         }
         else
         {
             StsException.systemError("StsEdgeLoopRadialGridLink.getLinkIndex() failed. Undefined arguments: " +
                 " rowOrCol: " + rowOrCol + " direction: " + direction);
             return NO_LINK;
         }
     }

     static public int getGridLinkIndex(int rowOrCol, int direction)
     {
         if (direction == NONE) return NO_LINK;

         if (rowOrCol == NONE)
         {
             if (direction == PLUS) return NEXT_POINT;
             else return PREV_POINT;
         }
         else if (rowOrCol == ROW)
         {
             if (direction == PLUS) return GRID_ROW_PLUS;
             else return GRID_ROW_MINUS;
         }
         else if (rowOrCol == COL)
         {
             if (direction == PLUS) return GRID_COL_PLUS;
             else return GRID_COL_MINUS;
         }
         else
         {
             StsException.systemError("StsEdgeLoopRadialGridLink.getLinkIndex() failed. Undefined arguments: " +
                 " rowOrCol: " + rowOrCol + " direction: " + direction);
             return NO_LINK;
         }
     }

     public byte getConnectType(int rowOrCol, int direction)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         if (linkIndex == NO_LINK) return NOT_CONNECTED;
         return connectTypes[linkIndex];
     }

     public byte getConnectType(int linkIndex)
     {
         return connectTypes[linkIndex];
     }

     public void setConnectType(int rowOrCol, int direction, byte connectType)
     {
         setConnectType(getLinkDirection(rowOrCol, direction), connectType);
     }

     public void setConnectType(int linkIndex, byte connectType)
     {
         if (connectTypes.length <= linkIndex) return;
         connectTypes[linkIndex] = connectType;
     }

     public StsEdgeLoopRadialGridLink getLink(int rowOrCol, int direction)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         return pointLinks[linkIndex];
     }

     public StsEdgeLoopRadialGridLink getLink(int linkIndex)
     {
         if (pointLinks == null || pointLinks.length <= linkIndex) return null;
         return pointLinks[linkIndex];
     }

     public StsEdgeLoopRadialGridLink getConnectedLink(int rowOrCol, int direction)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         return getConnectedLink(linkIndex);
     }

     public StsEdgeLoopRadialGridLink getConnectedLink(int linkIndex)
     {
         StsEdgeLoopRadialGridLink nextLink = pointLinks[linkIndex];
         if (nextLink == null) return null;
         if (!isConnected(linkIndex)) return null;
         return nextLink;
     }

     public int[] getOutsideAdjacentRowAndCol()
     {
         int row, col;
         int outsideDirection;

         int rowOrCol = rowCol.getRowOrCol();
         float rowF = rowCol.getRowF();
         float colF = rowCol.getColF();

         if (rowOrCol == NONE)
             return null;
         else if (rowOrCol == ROWCOL)
             return new int[]{(int) rowF, (int) colF};
         else
         {
             outsideDirection = getOutsideDirection(ROW);
             if (outsideDirection == MINUS)
                 col = StsMath.below(colF);
             else
                 col = StsMath.above(colF);

             outsideDirection = getOutsideDirection(COL);
             if (outsideDirection == MINUS)
                 row = StsMath.below(rowF);
             else
                 row = StsMath.above(rowF);
         }

         return new int[]{row, col};
     }

     public void initialSetLink(int rowOrCol, int direction, StsEdgeLoopRadialGridLink otherLink)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         pointLinks[linkIndex] = otherLink;
     }

     public void setLink(int rowOrCol, int direction, StsEdgeLoopRadialGridLink otherLink)
     {
         setLink(getLinkDirection(rowOrCol, direction), otherLink);
     }


     public void setLink(int linkIndex, StsEdgeLoopRadialGridLink otherLink)
     {
         if (pointLinks.length <= linkIndex) resizePointLinks(10);
         pointLinks[linkIndex] = otherLink;
     }

     private void resizePointLinks(int newSize)
     {
         pointLinks = (StsEdgeLoopRadialGridLink[]) StsMath.arraycopy(pointLinks, newSize, StsEdgeLoopRadialGridLink.class);
         byte[] newConnectTypes = new byte[newSize];
         System.arraycopy(connectTypes, 0, newConnectTypes, 0, connectTypes.length);
         connectTypes = newConnectTypes;
         linkedDirections = null;
     }

     /*
         public void setLink(int linkIndex, StsEdgeLoopRadialGridLink otherLink)
            {
                pointLinks[linkIndex] = otherLink;
                if(otherLink != null)
                {
                    recomputeLinkedDirections();
                    setConnectType(linkIndex, INSIDE_CONNECTED);
                }
            }
        */
     public boolean constructLinkedDirections()
     {
         return checkConstructLinkedDirections(true);
     }

     // recompute array of indices from prevPoint CCW to nextPoint including pointLinks in between.
     public boolean checkConstructLinkedDirections(boolean includeEdgeLinks)
     {
         byte[] newLinkedDirections = constructLinkedDirections(includeEdgeLinks);
         if(newLinkedDirections == null) return false;
         linkedDirections = newLinkedDirections;
         return true;
     }

     public byte[] constructLinkedDirections(boolean includeEdgeLinks)
     {
         boolean checkAngle;
         float prevAngle = 0.0f, nextAngle = 4.0f;

         if (pointLinks.length < 6) return null;

         ArrayList<OrderedDirection> angleOrderedDirections = new ArrayList<OrderedDirection>();

         StsEdgeLoopRadialGridLink prevLink = pointLinks[PREV_POINT];
         StsEdgeLoopRadialGridLink nextLink = pointLinks[NEXT_POINT];

         if (prevLink != null && nextLink != null)
         {
             checkAngle = true;
             prevAngle = getQuadAngle(prevLink);
             nextAngle = getQuadAngle(nextLink);

             if (nextAngle < prevAngle) nextAngle += 4.0f;

             if (includeEdgeLinks)
             {
                 addAngleOrderedDirection(angleOrderedDirections, PREV_POINT, prevAngle);
                 addAngleOrderedDirection(angleOrderedDirections, NEXT_POINT, nextAngle);
             }
         }
         else // accept all angles
             checkAngle = false;

         for (byte i = 0; i < 4; i++)
         {
             StsEdgeLoopRadialGridLink link = pointLinks[i];
             if (link != null)
             {
                 //                if(link == prevLink || link == nextLink) continue;
                 float angle = (float) i;
                 if (checkAngle)
                 {
                     if (angle < prevAngle) angle += 4.0f;
                     if (includeEdgeLinks && angle > prevAngle && angle < nextAngle)
                     {
                         addAngleOrderedDirection(angleOrderedDirections, i, angle);
                         connectTypes[i] = INSIDE_CONNECTED;
                     }

                     if (!includeEdgeLinks && angle >= prevAngle && angle <= nextAngle)
                     {
                         addAngleOrderedDirection(angleOrderedDirections, i, angle);
                         connectTypes[i] = INSIDE_CONNECTED;
                     }
                     else
                     {
                         //						StsException.systemError("StsEdgeLoopRadialGridLink.constructLinkedDirections() failed." +
                         //							"Couldn't add insideLink: " + link.getLabel());
                         //						return false;
                     }
                 }
                 else
                     addAngleOrderedDirection(angleOrderedDirections, i, angle);
             }
         }
         if (pointLinks.length == 10)
         {
             for (byte i = 6; i < 10; i++)
             {
                 StsEdgeLoopRadialGridLink link = pointLinks[i];
                 if (link != null)
                 {
                     if (link == prevLink || link == nextLink) continue;
                     float angle = getQuadAngle(link);
                     if (checkAngle)
                     {
                         if (angle < prevAngle) angle += 4.0f;
                         if (includeEdgeLinks && angle > prevAngle && angle < nextAngle)
                         {
                             addAngleOrderedDirection(angleOrderedDirections, i, angle);
                             connectTypes[i] = INSIDE_CONNECTED;
                         }

                         if (!includeEdgeLinks && angle >= prevAngle && angle <= nextAngle)
                         {
                             addAngleOrderedDirection(angleOrderedDirections, i, angle);
                             connectTypes[i] = INSIDE_CONNECTED;
                         }
                         else
                         {
                             // addAngleOrderedDirection(angleOrderedDirections, i, angle);
                             //StsException.systemError("StsEdgeLoopRadialGridLink.constructLinkedDirections() failed." +
                             //                         "Couldn't add insideLink: " + link.toString());
                             // return null;
                         }
                     }
                     else
                         addAngleOrderedDirection(angleOrderedDirections, i, angle);
                 }
             }
         }
         Collections.sort(angleOrderedDirections);
         // Object[] orderedDirections = angleOrderedDirections.values().toArray();
         int nOrderedDirections = angleOrderedDirections.size();
         byte[] tempLinkedDirections = new byte[nOrderedDirections];
         for (int i = 0; i < nOrderedDirections; i++)
             tempLinkedDirections[i] = angleOrderedDirections.get(i).direction;
         setLinkedDirections(tempLinkedDirections);
         return tempLinkedDirections;
     }

     private void addAngleOrderedDirection(ArrayList<OrderedDirection> orderedDirections, byte direction, float angle)
     {
         orderedDirections.add(new OrderedDirection(direction, angle));
     }

     public void setLinkedDirections(byte[] directions)
     {
         linkedDirections = directions;
         int nDirections = directions.length;
         int nInsideDirections = nDirections - 2;
         if(nInsideDirections == 0)
             insideDirection = NO_LINK;
         else if(nInsideDirections == 1)
             insideDirection = linkedDirections[1];
         else
             insideDirection = MULTI_INSIDE_LINKS;
     }

     public void setLinkedDirections(byte[] directions, int nLinks)
     {
         linkedDirections = new byte[nLinks];
         System.arraycopy(directions, 0, linkedDirections, 0, nLinks);
         int nInsideDirections = nLinks - 2;
         if(nInsideDirections == 0)
             insideDirection = NO_LINK;
         else if(nInsideDirections == 1)
             insideDirection = linkedDirections[1];
         else
             insideDirection = MULTI_INSIDE_LINKS;
     }

     /*
     synchronized public void setLinkedDirections(byte[] directions)
         {
             if(directions == null)
             {
                 this.linkedDirections = null;
                 return;
             }
             int nLinks = directions.length;
             linkedDirections = null;
             byte[] localLinkedDirections = new byte[nLinks];
             System.arraycopy(directions, 0, localLinkedDirections, 0, nLinks);
             linkedDirections = localLinkedDirections;
             if(nLinks != linkedDirections.length)
             {
                 System.out.println("Length is: " + linkedDirections.length + " should be: " + nLinks);
             }
         }
     */
     class OrderedDirection implements Comparable<OrderedDirection>
     {
         byte direction;
         float angle;

         OrderedDirection(byte direction, float angle)
         {
             this.direction = direction;
             this.angle = angle;
         }

         public int compareTo(OrderedDirection other)
         {
             if (angle < other.angle) return -1;
             if (angle > other.angle) return 1;
             return 0;
         }
     }

     private void addOrderedDirection(TreeMap orderedDirections, byte direction, float angle)
     {
         Byte dirObject = new Byte(direction);
         Float angleObject = new Float(angle);
         orderedDirections.put(angleObject, dirObject);
     }
     /*
         public boolean constructLinkedDirections()
            {
                boolean checkAngle;
                float prevAngle = 0.0f, nextAngle = 4.0f;

                if(pointLinks.length < 6) return false;

                TreeMap angleOrderedDirections = new TreeMap();

                StsEdgeLoopRadialGridLink prevLink = pointLinks[PREV_POINT];
                StsEdgeLoopRadialGridLink nextLink = pointLinks[NEXT_POINT];

                if(prevLink != null && nextLink != null)
                {
                    checkAngle = true;
                    prevAngle = getQuadAngle(prevLink);
                    nextAngle = getQuadAngle(nextLink);

                    if(nextAngle < prevAngle) nextAngle += 4.0f;

                    addOrderedDirection(angleOrderedDirections, (byte)PREV_POINT, prevAngle);
                    addOrderedDirection(angleOrderedDirections, (byte)NEXT_POINT, nextAngle);
                }
                else // accept all angles
                    checkAngle = false;

                for(byte i = 0; i < 4; i++)
                {
                    StsEdgeLoopRadialGridLink link = pointLinks[i];
                    if(link != null)
                    {
        //                if(link == prevLink || link == nextLink) continue;
                        float angle = (float)i;
                        if(checkAngle)
                        {
                            if(angle < prevAngle) angle += 4.0f;
                            if(angle < nextAngle)
                            {
                                addOrderedDirection(angleOrderedDirections, i, angle);
                                connectTypes[i] = INSIDE_CONNECTED;
                            }
                            else
                            {
        //						StsException.systemError("StsEdgeLoopRadialGridLink.constructLinkedDirections() failed." +
        //							"Couldn't add insideLink: " + link.getLabel());
        //						return false;
                            }
                        }
                        else
                            addOrderedDirection(angleOrderedDirections, i, angle);
                    }
                }
                if(pointLinks.length == 10)
                {
                    for(byte i = 6; i < 10; i++)
                    {
                        StsEdgeLoopRadialGridLink link = pointLinks[i];
                        if(link != null)
                        {
                            if(link == prevLink || link == nextLink) continue;
                            float angle = getQuadAngle(link);
                            if(checkAngle)
                            {
                                if(angle <= prevAngle) angle += 4.0f;
                                if(angle < nextAngle)
                                {
                                    addOrderedDirection(angleOrderedDirections, i, angle);
                                    connectTypes[i] = INSIDE_CONNECTED;
                                }
                                else
                                {
        //							StsException.systemError("StsEdgeLoopRadialGridLink.constructLinkedDirections() failed." +
        //								"Couldn't add insideLink: " + link.getLabel());
        //							return false;
                                }
                            }
                            else
                                addOrderedDirection(angleOrderedDirections, i, angle);
                        }
                    }
                }
                Object[] orderedDirections = angleOrderedDirections.values().toArray();
                int nOrderedDirections = orderedDirections.length;
                linkedDirections  = new byte[nOrderedDirections];
                for(int i = 0; i < nOrderedDirections; i++)
                    linkedDirections[i] = ((Byte)orderedDirections[i]).byteValue();

                return true;
            }
        */

     /*
         public void recomputeLinkedDirections()
            {
                if(pointLinks.length < 6) return;

                StsEdgeLoopRadialGridLink prevLink = pointLinks[PREV_POINT];
                StsEdgeLoopRadialGridLink nextLink = pointLinks[NEXT_POINT];

                float prevAngle = getQuadAngle(prevLink);
                float nextAngle = getQuadAngle(nextLink);

                int linkStartIndex = StsMath.above(prevAngle);
                int linkEndIndex = StsMath.below(nextAngle);
                if(nextAngle < prevAngle) linkEndIndex += 4;

                byte[] tempLinkedDirections = new byte[6];

                int nLinks = 0;
                tempLinkedDirections[nLinks++] = PREV_POINT;
                for(int i = linkStartIndex; i <= linkEndIndex; i++)
                {
                    int direction = i%4;
                    if(pointLinks[direction] != null)
                    {
                        tempLinkedDirections[nLinks++] = (byte)direction;
                        connectTypes[direction] = INSIDE_CONNECTED;
                    }
                }
                tempLinkedDirections[nLinks++] = NEXT_POINT;

                linkedDirections  = new byte[nLinks];
                System.arraycopy(tempLinkedDirections, 0, linkedDirections, 0, nLinks);
            }
        */

     public int[] getNextRowColLinkTowards(StsEdgeLoopRadialGridLink otherLink)
     {
         StsRowCol otherRowCol = otherLink.getRowCol();

         float rowF = rowCol.getRowF();
         float colF = rowCol.getColF();
         float otherRowF = otherRowCol.getRowF();
         float otherColF = otherRowCol.getColF();

         if (rowF == otherRowF)
         {
             if (otherColF - colF > 1.0f)
                 return new int[]{(int) rowF, StsMath.above(colF)};
             else if (colF - otherColF > 1.0f)
                 return new int[]{(int) rowF, StsMath.below(colF)};
             else
                 return null;
         }
         else if (colF == otherColF)
         {
             if (otherRowF - rowF > 1.0f)
                 return new int[]{StsMath.above(rowF), (int) colF};
             else if (rowF - otherRowF > 1.0f)
                 return new int[]{StsMath.below(rowF), (int) colF};
             else
                 return null;
         }
         else
             return null;
     }

     // check that links is used only twice if this is an edgeLink
     public boolean usageOK()
     {
         if (pointLinks.length < 6) return true;
         return ++flag <= 2;
     }

     public void debugDisplayLinkConnects(StsGLPanel3d glPanel3d)
     {
         StsColor color = StsColor.CYAN;
         float xyz[] = point.getPoint().v;
         StsColor connectColor;
         double leftRightShift, upDownShift;
         int size = 4;

         // draw center with rowOrCol color
         StsGLDraw.drawPoint(xyz, color, glPanel3d, size, size, 0.0);

         if (connectTypes[0] != NOT_CONNECTED && connectTypes[0] != DONT_CONNECT)
         {
             connectColor = getConnectColor(connectTypes[0]);
             leftRightShift = size;
             upDownShift = 0;
             StsGLDraw.drawPoint(xyz, connectColor, glPanel3d, size, size, 0.0, leftRightShift, upDownShift);
         }
         if (connectTypes[1] != NOT_CONNECTED && connectTypes[1] != DONT_CONNECT)
         {
             connectColor = getConnectColor(connectTypes[1]);
             leftRightShift = 0;
             upDownShift = size;
             StsGLDraw.drawPoint(xyz, connectColor, glPanel3d, size, size, 0.0, leftRightShift, upDownShift);
         }
         if (connectTypes[2] != NOT_CONNECTED && connectTypes[2] != DONT_CONNECT)
         {
             connectColor = getConnectColor(connectTypes[2]);
             leftRightShift = -size;
             upDownShift = 0;
             StsGLDraw.drawPoint(xyz, connectColor, glPanel3d, size, size, 0.0, leftRightShift, upDownShift);
         }
         if (connectTypes[3] != NOT_CONNECTED && connectTypes[3] != DONT_CONNECT)
         {
             connectColor = getConnectColor(connectTypes[3]);
             leftRightShift = 0;
             upDownShift = -size;
             StsGLDraw.drawPoint(xyz, connectColor, glPanel3d, size, size, 0.0, leftRightShift, upDownShift);
         }
     }

     public void debugDisplay(StsGLPanel3d glPanel3d, float offSet)
     {
         boolean minusConnect, plusConnect;
         StsColor stsColor;

         float[] pv = getXYZ();
         float[] v = new float[]{pv[0], pv[1], pv[2]};

         int rowOrCol = rowCol.getRowOrCol();

         if (rowOrCol == ROW)
         {
             stsColor = StsColor.CYAN;
             minusConnect = isConnected(ROW, MINUS);
             plusConnect = isConnected(ROW, PLUS);
             if (minusConnect && !plusConnect)
             {
                 if (offSet != 0.0f) v[0] -= offSet;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 6, 3, 0.0, 2.0, 0.0);
             }
             else if (plusConnect && !minusConnect)
             {
                 if (offSet != 0.0f) v[0] += offSet;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 6, 0.0, -2.0, 0.0);
             }
             else
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 3, 0.0, 0.0, 0.0);
         }
         else if (rowOrCol == COL)
         {
             stsColor = StsColor.MAGENTA;
             minusConnect = isConnected(COL, MINUS);
             plusConnect = isConnected(COL, PLUS);
             if (minusConnect && !plusConnect)
             {
                 if (offSet != 0.0f) v[1] -= offSet;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 6, 3, 0.0, 2.0, 0.0);
             }
             else if (plusConnect && !minusConnect)
             {
                 if (offSet != 0.0f) v[1] += offSet;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 6, 0.0, -2.0, 0.0);
             }
             else
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 3, 0.0, 0.0, 0.0);
         }
         else
         {
             stsColor = StsColor.WHITE;
             StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 3, 0.0, 0.0, 0.0);
         }
     }

     public void debugDisplay(StsGLPanel3d glPanel3d, StsSpectrum spectrum, boolean displayID, int seqNo, int rowOrCol)
     {
         boolean minusConnect, plusConnect;
         StsColor stsColor;
         float[] v = new float[3];

         float[] pv = getXYZ();

         stsColor = spectrum.getColor(seqNo % 8);

         v[2] = pv[2];

         if (rowOrCol == ROW)
         {
             minusConnect = isConnected(ROW, MINUS);
             plusConnect = isConnected(ROW, PLUS);
             v[1] = pv[1];
             if (minusConnect && !plusConnect)
             {
                 v[0] = pv[0] - 50.0f;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 6, 3, 0.0, 2.0, 0.0);
             }
             else if (plusConnect && !minusConnect)
             {
                 v[0] = pv[0] + 50.0f;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 6, 0.0, -2.0, 0.0);
             }
             else
             {
                 v[0] = pv[0];
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 3, 0.0, 0.0, 0.0);
             }
         }
         else if (rowOrCol == COL)
         {
             minusConnect = isConnected(COL, MINUS);
             plusConnect = isConnected(COL, PLUS);
             v[0] = pv[0];
             if (minusConnect && !plusConnect)
             {
                 v[1] = pv[1] - 50.0f;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 6, 3, 0.0, 2.0, 0.0);
             }
             else if (plusConnect && !minusConnect)
             {
                 v[1] = pv[1] + 50.0f;
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 6, 0.0, -2.0, 0.0);
             }
             else
             {
                 v[1] = pv[1];
                 StsGLDraw.drawPoint(v, stsColor, glPanel3d, 3, 3, 0.0, 0.0, 0.0);
             }
         }
     }

     private StsColor getConnectColor(byte connectType)
     {
         switch (connectType)
         {
             case NOT_CONNECTED:
                 return StsColor.GREY;
             case EDGE_CONNECTED:
                 return StsColor.MAGENTA;
             case INSIDE_CONNECTED:
                 return StsColor.YELLOW;
             case OUTSIDE_CONNECTED:
                 return StsColor.BLUE;
             case NEED_CONNECT:
                 return StsColor.WHITE;
             case DONT_CONNECT:
                 return StsColor.RED;
             default:
                 return StsColor.GREY;
         }
     }

     public int getNeedHasConnectDirection(int rowOrCol)
     {
         byte plusConnect = getConnectType(rowOrCol, PLUS);
         byte minusConnect = getConnectType(rowOrCol, MINUS);

         boolean needHasPlusConnect = (plusConnect != NOT_CONNECTED);
         boolean needHasMinusConnect = (minusConnect != NOT_CONNECTED);

         if (needHasPlusConnect && !needHasMinusConnect) return PLUS;
         else if (!needHasPlusConnect && needHasMinusConnect) return MINUS;
         else return NONE;
     }

     public int getEdgeDirection(int rowOrCol)
     {
         boolean plusConnected = isConnected(rowOrCol, PLUS);
         boolean minusConnected = isConnected(rowOrCol, MINUS);

         if (plusConnected && !minusConnected) return PLUS;
         else if (minusConnected && !plusConnected) return MINUS;
         else return NONE;
     }

     public StsEdgeLoopRadialGridLink getNextCounterClockwiseLink(StsEdgeLoopRadialGridLink prevLink)
     {
         if (pointLinks == null) return null;

         try
         {
             if (pointLinks.length == 4)
             {
                 for (int n = 0; n < 4; n++)
                 {
                     if (prevLink == pointLinks[n % 4])
                     {
                         for (int nn = n + 1; nn < n + 4; nn++)
                         {
                             StsEdgeLoopRadialGridLink nextLink = pointLinks[nn % 4];
                             if (nextLink != null) return nextLink;
                         }
                         return null;
                     }
                 }
                 return null;
             }

             if (!constructLinkedDirections()) return null;
             //			if(linkedDirections == null && !constructLinkedDirections()) return null;

             int nLinks = linkedDirections.length;
             for (int n = 0; n < nLinks; n++)
                 if (prevLink == pointLinks[linkedDirections[n]])
                     return pointLinks[linkedDirections[(n + 1) % nLinks]];

             return null;
         }
         catch (Exception e)
         {
             StsException.outputException("StsEdgeLoopRadialGridLink.getNextCounterClockwiseLink() failed.",
                 e, StsException.WARNING);
             return null;
         }
     }

     public StsEdgeLoopRadialGridLink getNextCounterClockwiseLink(int rowOrCol, int direction)
     {
         StsEdgeLoopRadialGridLink prevLink = getConnectedLink(rowOrCol, direction);
         return getNextCounterClockwiseLink(prevLink);
     }

     public StsEdgeLoopRadialGridLink getNextCounterClockwiseLink(int linkIndex)
     {
         StsEdgeLoopRadialGridLink prevLink = getConnectedLink(linkIndex);
         return getNextCounterClockwiseLink(prevLink);
     }

     /*
      public StsEdgeLoopRadialGridLink getNextCounterClockwisePolygonLink(StsEdgeLoopRadialGridLink prevLink)
         {
             int linkIndex = getCounterClockwiseLinkIndex(prevLink);
             if(linkIndex == NO_LINK) return null;
             if(connectTypes[linkIndex] != FAR_CONNECTED)
                 return pointLinks[linkIndex];
             else
                 return null;
         }
     */
     public int getCounterClockwiseLinkIndex(StsEdgeLoopRadialGridLink prevLink)
     {
         if (pointLinks == null) return NO_LINK;

         if (pointLinks.length == 4)
         {
             for (int n = 0; n < 4; n++)
             {
                 if (prevLink == pointLinks[n]) return (n + 1) % 4;
             }
             return NO_LINK;
         }
         else
         {
             byte[] linkedDirections = constructLinkedDirections(true);
             if (linkedDirections == null) return NO_LINK;
             //			if(linkedDirections == null && !constructLinkedDirections()) return NO_LINK;

             int nLinks = linkedDirections.length;
             for (int n = 0; n < nLinks; n++)
                 if (prevLink == pointLinks[linkedDirections[n]])
                     return linkedDirections[(n + 1) % nLinks];

             return NO_LINK;
         }
     }

     public StsEdgeLoopRadialGridLink getOppositeLink(StsEdgeLoopRadialGridLink prevLink)
     {
         for (int n = 0; n < 4; n++)
             if (prevLink == pointLinks[n])
                 return pointLinks[(n + 2) % 4];

         if (pointLinks.length < 10) return null;

         if (pointLinks[6] == prevLink) return pointLinks[8];
         if (pointLinks[7] == prevLink) return pointLinks[9];
         if (pointLinks[8] == prevLink) return pointLinks[6];
         if (pointLinks[9] == prevLink) return pointLinks[7];

         return null;
     }

     // a connection is ok if both types are DONT_CONNECT or both are not DONT_CONNECT
     protected boolean isConnectionOK(byte connectType1, byte connectType2)
     {
         if (connectType1 == DONT_CONNECT) return connectType2 == DONT_CONNECT;
         if (connectType2 == DONT_CONNECT) return connectType1 == DONT_CONNECT;
         return true;
     }

     public void connectToNextLink(int rowOrCol, StsEdgeLoopRadialGridLink nextLink, StsSurfaceGridable grid)
     {
         int plusLinkIndex = getLinkDirection(rowOrCol, PLUS);
         int minusLinkIndex = getLinkDirection(rowOrCol, MINUS);
         byte connectType = connectTypes[plusLinkIndex];
         byte nextConnectType = nextLink.getConnectType(minusLinkIndex);

         if (!isConnectionOK(connectType, nextConnectType))
         {
             StsException.systemError("StsEdgeLoopRadialGridLink.connectToNextLink() failed. " +
                 point.getLabel() + " connectType: " + connectType + " " +
                 nextLink.getPoint().getLabel() + " connectType: " + nextConnectType);
             return;
         }

         setLink(plusLinkIndex, nextLink);
         nextLink.setLink(minusLinkIndex, this);

         connectType = computeConnectType(rowOrCol, PLUS, grid);
         setConnectType(plusLinkIndex, connectType);
         nextLink.setConnectType(minusLinkIndex, connectType);
     }

     private byte computeConnectType(int rowOrCol, int direction, StsSurfaceGridable grid)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         if (linkIndex == NO_LINK) return NOT_CONNECTED;

         byte connectType = connectTypes[linkIndex];
         if (connectType == DONT_CONNECT) return DONT_CONNECT;

         StsEdgeLoopRadialGridLink otherLink = pointLinks[linkIndex];

         if (otherLink == null)
             return NOT_CONNECTED;

         return INSIDE_CONNECTED;
     }

     // Add a link from this to nextLink along this rowOrCol in the PLUS direction and
     // a MINUS connection from nextLink back to this link.
     public void connectToLinks(int rowOrCol, StsEdgeLoopRadialGridLink prevLink, StsEdgeLoopRadialGridLink nextLink, StsSurfaceGridable grid)
     {
         int plusLinkIndex = getLinkDirection(rowOrCol, PLUS);
         int minusLinkIndex = getLinkDirection(rowOrCol, MINUS);
         byte prevConnectType = prevLink.getConnectType(plusLinkIndex);
         byte nextConnectType = nextLink.getConnectType(minusLinkIndex);

         if (!isConnectionOK(prevConnectType, nextConnectType))
         {
             StsException.systemError("StsEdgeLoopRadialGridLink.connectToLinks() failed. " +
                 prevLink.toString() + " connectType: " + prevConnectType + " " +
                 nextLink.toString() + " connectType: " + nextConnectType);
             return;
         }

         byte connectType;

         prevLink.setLink(plusLinkIndex, this);
         setLink(minusLinkIndex, prevLink);
         connectType = computeConnectType(rowOrCol, MINUS, grid);
         prevLink.setConnectType(plusLinkIndex, connectType);
         setConnectType(minusLinkIndex, connectType);

         setLink(plusLinkIndex, nextLink);
         nextLink.setLink(minusLinkIndex, this);
         connectType = computeConnectType(rowOrCol, PLUS, grid);
         setConnectType(plusLinkIndex, connectType);
         nextLink.setConnectType(minusLinkIndex, connectType);
     }

     public StsEdgeLoopRadialGridLink getNextOutsideLink(int linkIndex)
     {
         StsEdgeLoopRadialGridLink link = null;
         StsEdgeLoopRadialGridLink nextLink = getConnectedLink(linkIndex);
         while (nextLink != null)
         {
             link = nextLink;
             nextLink = link.getConnectedLink(linkIndex);
         }
         return link;
     }

     public int[] getLowerLeftRowCol(StsEdgeLoopRadialGridLink nextLink, byte polygonType, StsSurfaceGridable grid)
     {
         StsRowCol rowCol = point.getRowCol(grid, polygonType);
         if (rowCol == null)
         {
             StsException.systemError(this, "getLowerLeftGridRowCol", "StsGridRowCol is null for link " + toString());
             return new int[]{0, 0};
         }
         StsRowCol nextRowCol = nextLink.point.getRowCol(grid, polygonType);
         return rowCol.getLowerLeftRowCol(nextRowCol);
     }

     public void delete()
     {
         StsEdgeLoopRadialGridLink link;
         int nLinks = pointLinks.length;
         for (int n = 0; n < nLinks; n++)
         {
             link = pointLinks[n];
             if (link != null)
             {
                 removeLink(link);
                 link.removeLink(this);
             }
         }
     }

     protected void removeLink(StsEdgeLoopRadialGridLink link)
     {
         int nLinks = pointLinks.length;
         for (int n = 0; n < nLinks; n++)
         {
             if (pointLinks[n] == link)
             {
                 pointLinks[n] = null;
                 connectTypes[n] = NOT_CONNECTED;
             }
         }
     }

     // return the direction from this fault point onto the grid
     public int getFaultConnectDirection(int rowOrCol)
     {
         if (rowOrCol == ROW)
         {
             if (connectTypes[ROW_PLUS] == DONT_CONNECT &&
                 pointLinks[ROW_MINUS] != null) return ROW_MINUS;
             else if (connectTypes[ROW_MINUS] == DONT_CONNECT &&
                 pointLinks[ROW_PLUS] != null) return ROW_PLUS;
         }
         else if (rowOrCol == COL)
         {
             if (connectTypes[COL_PLUS] == DONT_CONNECT &&
                 pointLinks[COL_MINUS] != null) return COL_MINUS;
             else if (connectTypes[COL_MINUS] == DONT_CONNECT &&
                 pointLinks[COL_PLUS] != null) return COL_PLUS;
         }
         return NONE;
     }

     public int getFaultDirection(int rowOrCol)
     {
         if (getConnectType(rowOrCol, PLUS) == DONT_CONNECT)
             return MINUS;
         else if (getConnectType(rowOrCol, MINUS) == DONT_CONNECT)
             return PLUS;
         else
             return NONE;
     }

     public StsEdgeLoopRadialGridLink getLastConnectedLink(int rowOrCol, int direction)
     {
         return getLastConnectedLink(getLinkDirection(rowOrCol, direction));
     }

     // Search in this direction for last connected link
     public StsEdgeLoopRadialGridLink getLastConnectedLink(int linkIndex)
     {
         StsEdgeLoopRadialGridLink nextLink = getConnectedLink(linkIndex);
         if (nextLink == null) return null;

         StsEdgeLoopRadialGridLink link = nextLink;
         while ((nextLink = link.getConnectedLink(linkIndex)) != null)
             link = nextLink;
         return link;
     }

     public StsEdgeLoopRadialGridLink getLastInsideConnectedLink(int rowOrCol, int direction)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         StsEdgeLoopRadialGridLink nextLink = getInsideIndexedGridLink(linkIndex);
         if (nextLink == null) return null;
         if (!isConnected(linkIndex)) return null;
         StsEdgeLoopRadialGridLink link = nextLink;
         while ((nextLink = link.getConnectedLink(linkIndex)) != null)
             link = nextLink;
         return link;
     }

     public StsEdgeLoopRadialGridLink getNextInsideConnectedLink(int rowOrCol, int direction)
     {
         int linkIndex = getLinkDirection(rowOrCol, direction);
         StsEdgeLoopRadialGridLink nextLink = getInsideIndexedGridLink(linkIndex);
         if (nextLink == null) return null;
         if (!isConnected(linkIndex)) return null;
         return nextLink;
     }

     public StsEdgeLoopRadialGridLink[] getBoundingEdgeLinks(int rowOrCol)
     {
         StsEdgeLoopRadialGridLink link;

         if (!isRowOrCol()) return null;
         StsEdgeLoopRadialGridLink[] boundingEdgeLinks = new StsEdgeLoopRadialGridLink[2];
         boundingEdgeLinks[0] = getLink(rowOrCol, MINUS);
         boundingEdgeLinks[1] = getLink(rowOrCol, PLUS);
         return boundingEdgeLinks;
     }
     public StsSectionRowCol getSectionRowCol(StsSection section)
     {
         if(point == null) return null;
         return (StsSectionRowCol)point.getRowCol(section);
     }

     static public String getLabel(StsEdgeLoopRadialGridLink link)
     {
         if (link != null) return link.toString();
         else return new String(" is null. ");
     }

     public String toString()
     {
         int rowOrCol = rowCol.getRowOrCol();
         return new String("StsEdgeLoopRadialGridLink: " + StsParameters.rowCol(rowOrCol) +
             " rowF: " + rowCol.getRowF() + " colF: " + rowCol.getColF() +
             " inside direction: " + StsParameters.getGridLinkName(insideDirection));
     }

     static public void main(String[] args)
     {
        StsEdgeLoopRadialGridLink link = new StsEdgeLoopRadialGridLink();
            link.linkedDirections = new byte[0];
         byte[] directions = new byte[] { 1, 2, 0, 0, 0, 0 };
         link.setLinkedDirections(directions, 2);
     }
}
