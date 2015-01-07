package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public interface IPolygon
{
    public String getLabel();
    public double[][] getPoints();
    public int getNPoints();
    public int getID();
    public void drawConcaveFailed(int error);
}