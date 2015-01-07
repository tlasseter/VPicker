package com.Sts.Types;

import com.Sts.Interfaces.*;
 import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Apr 23, 2010
  * Time: 10:14:57 AM
  * To change this template use File | Settings | File Templates.
  */

/**
 * Inner class defining useful information about points on grid.  Defined static so we can construct these
  * indenpendent of the enclosing class.
  */

 public class StsGridCrossingPoint
 {
     public StsPoint point;
     /** xyz and f values at this point */
     public float iF, jF;
     /** float grid index coordinates */
     public int rowOrCol;
     /** indicates point is on ROW I or COL J */
     public int rowConnect = NONE;
     /** indicates direction on ROW or COL away from RIGHT side of line */
     public int colConnect = NONE;
     /** indicates direction on ROW or COL away from RIGStsXYGridable e */
     public float gridLength;

    protected static final int MINUS = StsParameters.MINUS;
    protected static final int PLUS = StsParameters.PLUS;
    protected static final int NONE = StsParameters.NONE;

    protected static final int ROW = StsParameters.ROW;
    protected static final int COL = StsParameters.COL;
    protected static final int ROWCOL = StsParameters.ROWCOL;
     /** length in grid units from start of edge */

     public StsGridCrossingPoint(StsXYGridable grid, StsPoint point)
     {
         this.point = point;
         iF = grid.getRowCoor(point.v);
         jF = grid.getColCoor(point.v);
     }

     public StsGridCrossingPoint(StsGridPoint stsGridPoint)
     {
         point = new StsPoint(stsGridPoint.getPoint());
         iF = stsGridPoint.rowF;
         jF = stsGridPoint.colF;
         rowOrCol = stsGridPoint.rowOrCol;
     }

     public StsGridCrossingPoint(StsGridCrossingPoint gridPoint)
     {
         point = new StsPoint(gridPoint.point);
         iF = gridPoint.iF;
         jF = gridPoint.jF;
         rowOrCol = gridPoint.rowOrCol;
     }

     public StsGridCrossingPoint(float iF_, float jF_, StsXYGridable grid)
     {
         iF = iF_;
         jF = jF_;
         float x = grid.getXCoor(iF, jF);
         float y = grid.getYCoor(iF, jF);
         point = new StsPoint(x, y, 0.0f);
         setRowOrCol();
     }

     public StsGridCrossingPoint(float rowF, float colF, float[] xyz, StsXYGridable grid)
     {
         iF = rowF;
         jF = colF;
         float x = grid.getXCoor(iF, jF);
         float y = grid.getYCoor(iF, jF);
         point = new StsPoint(xyz);
         setRowOrCol();
     }

     public StsGridCrossingPoint(float rowColF, int rowOrCol, StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1, float ff, int connect)
     {
         this.rowOrCol = rowOrCol;

         if (rowOrCol == ROW)
         {
             iF = rowColF;
             jF = gridPoint0.jF + ff * (gridPoint1.jF - gridPoint0.jF);
             rowConnect = connect;
         }
         else
         {
             jF = rowColF;
             iF = gridPoint0.iF + ff * (gridPoint1.iF - gridPoint0.iF);
             colConnect = connect;
         }

         point = StsPoint.staticInterpolatePoints(gridPoint0.point, gridPoint1.point, ff);
         gridLength = gridPoint0.gridLength + ff * (gridPoint1.gridLength - gridPoint0.gridLength);
     }

     public StsGridCrossingPoint(StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1, float ff)
     {
         rowOrCol = NONE;

         iF = gridPoint0.iF + ff * (gridPoint1.iF - gridPoint0.iF);
         jF = gridPoint0.jF + ff * (gridPoint1.jF - gridPoint0.jF);

         point = StsPoint.staticInterpolatePoints(gridPoint0.point, gridPoint1.point, ff);
         float z = point.getZorT();
         point.setZorT(StsMath.checkRoundOffInteger(z));
         gridLength = gridPoint0.gridLength + ff * (gridPoint1.gridLength - gridPoint0.gridLength);

         setRowOrCol();

         if (rowOrCol == ROW || rowOrCol == ROWCOL)
             rowConnect = StsMath.sign(gridPoint1.iF - gridPoint0.iF);
         if (rowOrCol == COL || rowOrCol == ROWCOL)
             colConnect = StsMath.sign(gridPoint0.jF - gridPoint1.iF);
     }

     public void computeGridLength(StsGridCrossingPoint otherPoint)
     {
         float dRow = iF - otherPoint.iF;
         float dCol = jF - otherPoint.jF;
         float dGridLength = (float) Math.sqrt(dRow * dRow + dCol * dCol);
         gridLength = otherPoint.gridLength + dGridLength;
     }

     void setRowOrCol()
     {
         boolean isRow = StsMath.isExactIntegral(iF);
         boolean isCol = StsMath.isExactIntegral(jF);

         //            if(isRow) iF = StsMath.roundOffInteger(iF);
         //            if(isCol) jF = StsMath.roundOffInteger(jF);

         if (isRow)
         {
             if (isCol) rowOrCol = ROWCOL;
             else rowOrCol = ROW;
         }
         else if (isCol)
             rowOrCol = COL;
         else
             rowOrCol = NONE;
     }

     void setConnects(int rowConnect, int colConnect)
     {
         if (rowOrCol == ROW || rowOrCol == ROWCOL)
             this.rowConnect = rowConnect;
         if (rowOrCol == COL || rowOrCol == ROWCOL)
             this.colConnect = colConnect;
     }

     public float getSectionF()
     {
         if (point.v.length < 4) return 0.0f;
         else return point.v[3];
     }

     public int[] getLowerLeftRowCol()
     {
         int row = (int) iF;
         int col = (int) jF;

         if (rowOrCol == ROW || rowOrCol == ROWCOL)
             if (rowConnect == MINUS) row--;
             else if (rowOrCol == COL || rowOrCol == ROWCOL)
                 if (colConnect == MINUS) col--;
         return new int[]{row, col};
     }

     public void adjustToCellGridding()
     {
         iF += 0.5f;
         jF += 0.5f;
         setRowOrCol();
     }

     public boolean sameAs(StsGridCrossingPoint otherGridPoint)
     {
         if (rowOrCol != otherGridPoint.rowOrCol) return false;
         if (rowOrCol == ROW || rowOrCol == ROWCOL)
             if ((int) iF != (int) otherGridPoint.iF) return false;
         if (rowOrCol == COL || rowOrCol == ROWCOL)
             if ((int) jF != (int) otherGridPoint.jF) return false;
         return true;
     }

     public String toString()
     {
         return "RowOrCol " + StsParameters.getRowOrColString(rowOrCol) +  " rowF " + iF + " col" + jF + " point " + point.toString();
     }
 }
