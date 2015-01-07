package com.Sts.DBTypes;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 5:51:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFaultStickSetClass extends StsClass
{
    int nextColorIndex = 0;

    public StsFaultStickSetClass()
    {
        userName = "Set of Fault Sticks";                                                              
    }

    public StsColor getNextColor()
    {
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        return new StsColor(spectrum.getColor(nextColorIndex++));
    }
}
