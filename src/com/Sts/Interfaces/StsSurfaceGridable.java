package com.Sts.Interfaces;

import com.Sts.Types.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public interface StsSurfaceGridable extends StsGridable
{
    float getRowCoor(float[] xyz);
    float getColCoor(float[] xyz);

    public StsPoint getPoint(int row, int col);
    public StsPoint getPoint(float rowF, float colF);
    public float[] getXYZorT(int row, int col);
    public float[] getXYZorT(float rowF, float colF);
    public float[] getNormal(int row, int col);
    float[] getNormal(float rowF, float colF);
    public void checkConstructGridNormals();
    public String toString();

}
