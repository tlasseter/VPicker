//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.UI.Beans.*;

abstract public class StsMainTimeObject extends StsMainObject
{

    public long bornDate = -1L;
    public long deathDate = -1l;

    public StsMainTimeObject()
    {
        super();
    }
    public StsMainTimeObject(boolean persistent)
    {
        super(persistent);
    }
    public StsMainTimeObject(boolean persistent, String name)
    {
        super(persistent);
        setName(name);
    }
    public boolean isAlive(long time)
    {
        if(bornDate <= 0)      // Time not set
            return true;

        // Time check
        if((bornDate < time) && ((deathDate > time) || (deathDate <= 0)))
            return true;
        return false;
    }
    public long getBornDateLong()
    {
        return bornDate;
    }
    public long getDeathDateLong()
    {
        return deathDate;
    }
    public String getBornDate()
    {
        if(bornDate < 0)
            return "Undefined";
        return StsDateFieldBean.convertToString(bornDate);
    }

    public void setBornDate(String born)
    {
        long newDate = 0l;
        try { newDate = StsDateFieldBean.convertToLong(born); }
        catch(Exception ex)  {  }
        setBornDate(newDate);
        return;
    }
    
    public void setDeathDate(long death)
    {
        if(death >= -1l)
        {
            deathDate = death;
            // deathDateString = StsDateFieldBean.convertToString(death);
            dbFieldChanged("deathDate", deathDate);            
        }
        return;
    }
    
    public void setBornDate(long born)
    {
        if(born >= -1l)
        {
            bornDate = born;
            // bornDateString = StsDateFieldBean.convertToString(born);
            dbFieldChanged("bornDate", bornDate);
        }
        return;
    }
    public String getDeathDate()
    {
        if(deathDate < 0)
            return "Undefined";
        return StsDateFieldBean.convertToString(deathDate);
    }

    public void setDeathDate(String death)
    {
        long newDate = 0l;
        try { newDate = StsDateFieldBean.convertToLong(death); }
        catch(Exception ex) { }
        setDeathDate(newDate);
        return;
    }


}
