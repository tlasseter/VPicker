package com.Sts.Actions.Wizards.Well;

import com.Sts.DBTypes.*;

public class StsWellFactory extends StsAbstractWellFactory
{
    static StsWellFactory wellFactory = null;
    static boolean isLive = false;

    public StsWellFactory()
    {
    }

    static StsWellFactory getWellFactory()
    {
        if(wellFactory == null) wellFactory = new StsWellFactory();
        return wellFactory;
    }

    /**
     * createWellInstance
     *
     * @return StsWell
     * @todo Implement this com.Sts.DBTypes.StsAbstractWellFactory method
     */
    public StsWell createWellInstance()
    {
        if(isLive)
            return new StsLiveWell();
        else
            return new StsWell();
    }
    public StsWell createWellInstance(String wellname)
    {
        if(isLive)
            return new StsLiveWell(wellname, true);
        else
            return new StsWell(wellname, true);
    }
    public void setIsLive(boolean live)
    {
        isLive = live;
    }
}
