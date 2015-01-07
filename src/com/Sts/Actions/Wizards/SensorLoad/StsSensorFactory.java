package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.DBTypes.*;

public class StsSensorFactory extends StsAbstractSensorFactory
{
    static StsSensorFactory sensorFactory = null;

    public StsSensorFactory()
    {
    }

    static StsSensorFactory getSensorFactory()
    {
        if(sensorFactory == null) sensorFactory = new StsSensorFactory();
        return sensorFactory;
    }

    /**
     * createSensorInstance
     *
     * @return StsSensor
     * @todo Implement this com.Sts.DBTypes.StsAbstractSensorFactory method
     */
    public StsSensor createSensorInstance()
    {
        return new StsSensor();
    }

    public StsSensor createSensorInstance(String sensorName, byte type)
    {
        if(type == StsSensor.DYNAMIC)
            return new StsDynamicSensor(null, sensorName);
        else
            return new StsStaticSensor(null, sensorName);
    }
}
