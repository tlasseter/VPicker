package com.Sts.Types.Simulation;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 7, 2010
 * Time: 3:40:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsNNC implements Comparable<StsNNC>
{
    public StsBlock.BlockCellColumn.GridCell rightGridCell;
    public StsBlock.BlockCellColumn.GridCell leftGridCell;
    private double[] faceCenter;
    private double[] faceNormal;
    /** initially this is the intersection area of the two polygons; subsequently converted to area normal to line between cell centers */
    public double area;
    public double areaFactor;
    public double trans;
    public double kEff;
    double lineLength;

    public StsNNC parentNNC = null;

    static final boolean debug = true;
    static final int[] debugIJK = new int[] { 66, 0, 0};
    static boolean debugIJKB;

    static public StsNNC constructor(float[] faceCenterF, double[] faceNormal, StsBlock.BlockCellColumn.GridCell rightGridCell,
                                     StsBlock.BlockCellColumn.GridCell leftGridCell, double area)
    {
        try
        {
            if(rightGridCell.volume == 0.0f)
            {
                StsException.systemDebug(StsNNC.class, "constructor", "bad NNC: cell volume is 0.0" + rightGridCell.toString());
                return null;
            }
            if(leftGridCell.volume == 0.0f)
            {
                StsException.systemDebug(StsNNC.class, "constructor", "bad NNC: cell volume is 0.0" + leftGridCell.toString());
                return null;
            }
            StsNNC nnc = new StsNNC(faceCenterF, faceNormal, rightGridCell, leftGridCell, area);
            // nnc.computeTrans();
            rightGridCell.addNNC(nnc);
            leftGridCell.addNNC(nnc);
            return nnc;
        }
        catch (StsException e)
        {
            StsException.systemError(StsNNC.class, "constructor", "Failed NNC connection between " + rightGridCell.toString() + " and " + leftGridCell.toString());
            return null;
        }
    }


    private StsNNC(float[] faceCenterF, double[] faceNormal, StsBlock.BlockCellColumn.GridCell rightGridCell, StsBlock.BlockCellColumn.GridCell leftGridCell, double area)
        throws StsException
    {
        this(StsMath.convertFloatToDoubleArray(faceCenterF), faceNormal, rightGridCell, leftGridCell, area);
    }

    private StsNNC(double[] faceCenter, double[] faceNormal, StsBlock.BlockCellColumn.GridCell rightGridCell, StsBlock.BlockCellColumn.GridCell leftGridCell, double area)
        throws StsException
    {
        if(StsBlock.debugIJK(leftGridCell))
            StsException.systemDebug(StsNNC.class, "constructor", "GridCell: " + leftGridCell.toIJKBString());
        if(StsBlock.debugIJK(rightGridCell))
            StsException.systemDebug(StsNNC.class, "constructor", "GridCell: " + rightGridCell.toIJKBString());
        this.rightGridCell = rightGridCell;
        this.leftGridCell = leftGridCell;
        this.faceCenter = faceCenter;
        this.faceNormal = faceNormal;
        this.area = area;
        if(area < 0.0)
            StsException.systemError(this, "constructor", "area < 0.0");
        if(Double.isNaN(area))
            StsException.systemError(this, "constructor", "area is NaN");
    }

    public StsNNC constructParentNNC()
    {
        try
        {
            StsBlock.BlockCellColumn.GridCell rightParentCell = rightGridCell.parentCell;
            StsBlock.BlockCellColumn.GridCell leftParentCell = leftGridCell.parentCell;
            if(rightParentCell == null && leftParentCell == null) return null;
            if(rightParentCell == null) rightParentCell = rightGridCell;
            if(leftParentCell == null) leftParentCell = leftGridCell;
            parentNNC = new StsNNC(faceCenter, faceNormal, rightParentCell, leftParentCell, area);
            // parentNNC.computeParentTrans();
            return parentNNC;
        }
        catch(Exception e)
        {
            StsException.systemError(this, "constructParentNNC", "Failed to construct parent NNC " + toString());
            return null;
        }
    }

    public void computeTrans() throws StsException
    {
        try
        {
            setSmallerIJKB();

            if(isDebug(this))
                StsException.systemDebug(toString());

            // compute intersection of line from right to left grid cell and polygon on fault section
            double[] intersection;

            double[] rightCenter = StsMath.convertFloatToDoubleArray(rightGridCell.cellCenter);
            double[] leftCenter = StsMath.convertFloatToDoubleArray(leftGridCell.cellCenter);
            double[] dLine = StsMath.subtract(leftCenter, rightCenter);

            lineLength = StsMath.length(dLine);
            if (lineLength <= 0.0)
            {
                trans = 0.0;
                return;
            }

            double distSqRight = StsMath.distanceSq(faceCenter,  rightCenter, 3);
            double distSqLeft = StsMath.distanceSq(faceCenter,  leftCenter, 3);
            double distSq = StsMath.distanceSq(rightCenter, leftCenter, 3);
            double fLeft = (1.0 + (distSqLeft - distSqRight)/distSq)/2;

            float[] rightK = rightGridCell.getKxyz();
            float[] leftK = leftGridCell.getKxyz();
            double kLeft = effectivePerm(leftK, dLine);
            double kRight = effectivePerm(rightK, dLine);
            double kWeight =  kLeft * (1 - fLeft) + kRight * fLeft;
            if (kWeight <= 0.0)
            {
                trans = 0.0;
                return;
            }
            kEff = kLeft * kRight / kWeight;
            // compute area faceNormal to line connecting centers
            areaFactor = Math.abs(StsMath.dot(dLine, faceNormal));
            //if(adjustArea)
                // this.area = area*areaFactor;
            trans = this.area * kEff / lineLength;
            if(trans < 0.0)
            {
                StsException.systemError(this, "computeTrans", "trans < 0.0");
                trans = 0.0;
            }
            if(Double.isNaN(trans) || Double.isInfinite(trans))
            {
                StsException.systemError(this, "constructor", "trans is Nan or Infinite at " + toString());
                trans = 0.0;
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeTrans", e);
            throw new StsException(StsException.WARNING, "NNC constructor failed.");
        }
    }

    public void adjustArea()
    {
        trans *= areaFactor;
        area *= areaFactor;
    }

    public StsNNC(StsBlock.BlockCellColumn.GridCell rightGridCell, StsBlock.BlockCellColumn.GridCell leftGridCell, double trans)
    {
        this.rightGridCell = rightGridCell;
        this.leftGridCell = leftGridCell;
        if(trans < 0.0)
        {
            StsException.systemError(this, "constructor", "trans < 0.0");
            trans = 0.0;
        }
        this.trans = trans;
        setSmallerIJKB();
    }

    static public boolean isDebug(StsNNC nnc)
    {
        return debug && (sameIJK(nnc.getIjkbLeft(), debugIJK) || sameIJK(nnc.getIjkbRight(), debugIJK));
    }

    private double effectivePerm(float[] perm, double[] vector)
    {
        double sum = 0.0;
        for (int n = 0; n < 3; n++)
        {
            if(perm[n] > 0.0)
                sum += vector[n] * vector[n] / perm[n];
        }
        if (sum <= 0.0) return 0;
        return 1 / sum;
    }

    public String toString()
    {
        return "right ijkb: " + rightGridCell.toIJKBString() + " " + "left ijkb: " + leftGridCell.toIJKBString() + " area: " + area + " TRANS: " + trans;
    }

    public String detailString()
    {
        return  " Keff: " + kEff + " length: " + lineLength + " areaFactor: " + areaFactor;
    }

    static public String getStringIJKB(int[] ijkb)
    {
        return ijkb[0] + " " + ijkb[1] + " " + ijkb[2] + " " + ijkb[3];
    }

    static public boolean ijkCompare(int[] ijkbA, int[] ijkbB)
    {
        for(int n = 0; n < 3; n++)
            if(ijkbA[n] != ijkbB[n]) return false;
        return true;
    }

    /** set smaller of two NNC indexSets to the "right" one to order for visual inspection. */
    void setSmallerIJKB()
    {
        for(int n = 3; n >= 0; n--)
        {
            if (getIjkbRight()[n] == getIjkbLeft()[n]) continue;
            if (getIjkbRight()[n] > getIjkbLeft()[n])
            {
                StsBlock.BlockCellColumn.GridCell tempGridCell = leftGridCell;
                leftGridCell = rightGridCell;
                rightGridCell = tempGridCell;
            }
            return;
        }
        StsException.systemError(this, "setSmallerIJKB", " NNC connects same two cells " + toString());
    }

    public int compareTo(StsNNC otherNNC)
    {
        for (int n = 3; n >= 0; n--)
        {
            int m = getIjkbRight()[n];
            int mm = otherNNC.getIjkbRight()[n];
            if (m < mm) return -1;
            else if (m > mm) return 1;
        }
        // both are the same, sort on IjkbLeft
        for (int n = 3; n >= 0; n--)
        {
            int m = getIjkbLeft()[n];
            int mm = otherNNC.getIjkbLeft()[n];
            if (m < mm) return -1;
            else if (m > mm) return 1;
        }
        return 0;
    }

    public boolean sameNNC(StsNNC otherNNC)
    {
        if(!sameIJKB(getIjkbLeft(), otherNNC.getIjkbLeft())) return false;
        if(!sameIJKB(getIjkbRight(), otherNNC.getIjkbRight())) return false;
        return true;
    }

    static  public boolean sameIJKB(int[] ijkbA, int[] ijkbB)
    {
        for (int n = 3; n >= 0; n--)
            if(ijkbA[n] != ijkbB[n]) return false;
        return true;
    }

    static public boolean sameIJK(int[] ijkbA, int[] ijkbB)
    {
        for (int n = 2; n >= 0; n--)
            if(ijkbA[n] != ijkbB[n]) return false;
        return true;
    }

    public String getEclipseRightBlockIndexString()
    {
        return getEclipseBlockIndexString(rightGridCell, getIjkbRight());
    }

    public String getEclipseLeftBlockIndexString()
    {
        return getEclipseBlockIndexString(leftGridCell, getIjkbLeft());
    }

    String getEclipseBlockIndexString(StsBlock.BlockCellColumn.GridCell gridCell, int[] ijkb)
    {
        StsBlock block = gridCell.getBlock();
        int row = block.eclipseRowMin + block.rowMax - ijkb[0];
        int col = ijkb[1] - block.colMin + 1;
        int layer = ijkb[2] + 1;
        return col + "  " + row + "  " + layer;
    }

    public int[] getIjkbRight()
    {
        return rightGridCell.getIJKB();
    }

    public int[] getIjkbLeft()
    {
        return leftGridCell.getIJKB();
    }
}