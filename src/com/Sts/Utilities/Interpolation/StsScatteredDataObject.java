package com.Sts.Utilities.Interpolation;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 24, 2007
 * Time: 11:57:43 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsScatteredDataObject
{
    public double distSq;
    public double weight;
    public Object valueObject = null;

    abstract public StsScatteredDataObject interpolateValue(StsScatteredDataObject result);

    public StsScatteredDataObject(double distSq, Object value)
    {
        this.distSq = distSq;
        this.valueObject = value;
    }

    public int compareTo(Object other)
    {
        StsScatteredDataObject otherDataObject = (StsScatteredDataObject)other;
        if(distSq > otherDataObject.distSq) return 1;
        else if(distSq < otherDataObject.distSq) return -1;
        else return 0;
    }

    public double getDistance() { return Math.sqrt(distSq); }
}
