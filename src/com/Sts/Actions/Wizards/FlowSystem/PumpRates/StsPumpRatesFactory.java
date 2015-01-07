package com.Sts.Actions.Wizards.FlowSystem.PumpRates;

import com.Sts.DBTypes.*;

public class StsPumpRatesFactory extends StsAbstractSensorFactory
{
    static StsPumpRatesFactory pumpFactory = null;

    public StsPumpRatesFactory()
    {
    }

    static StsPumpRatesFactory getTankFactory()
    {
        if(pumpFactory == null) pumpFactory = new StsPumpRatesFactory();
        return pumpFactory;
    }

    /**
     * createSensorInstance
     *
     * @return StsSensor
     * @todo Implement this com.Sts.DBTypes.StsAbstractSensorFactory method
     */
    public StsPump createSensorInstance()
    {
        return new StsPump();
    }
    public StsPump createSensorInstance(String pumpName, byte type)
    {
        return new StsPump(null, pumpName);
    }
}
