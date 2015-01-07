package com.Sts.Interfaces;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public interface StsVolumeDisplayable
{
    public byte[] readBytePlaneData(int dir, float dirCoordinate);
    public float getScaledValue(byte byteValue);
    public boolean isByteValueNull(byte byteValue);
}
