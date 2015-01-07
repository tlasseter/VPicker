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


public abstract class StsAbstractWellFactory
{
    abstract public StsWell createWellInstance();
    abstract public StsWell createWellInstance(String wellname);
    abstract public void setIsLive(boolean live);    
}
