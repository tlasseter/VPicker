package com.Sts.DBTypes;

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


public abstract class StsAbstractSensorFactory
{
    abstract public StsSensor createSensorInstance();
    abstract public StsSensor createSensorInstance(String sensorName, byte type);
}
