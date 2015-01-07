package com.Sts.Actions.Wizards.FlowSystem.TankLevels;

import com.Sts.DBTypes.*;

public class StsTankLevelsFactory extends StsAbstractSensorFactory
{
    static StsTankLevelsFactory tankFactory = null;

    public StsTankLevelsFactory()
    {
    }

    static StsTankLevelsFactory getTankFactory()
    {
        if(tankFactory == null) tankFactory = new StsTankLevelsFactory();
        return tankFactory;
    }

    /**
     * createSensorInstance
     *
     * @return StsSensor
     * @todo Implement this com.Sts.DBTypes.StsAbstractSensorFactory method
     */
    public StsTank createSensorInstance()
    {
        return new StsTank();
    }
    public StsTank createSensorInstance(String tankName, byte type)
    {
        return new StsTank(null, tankName, true);
    }
}
