package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.UI.Beans.*;

public class StsMovieClass extends StsClass implements StsSerializable
{

    boolean defaultLoop = true;
    boolean defaultCycleVolumes = true;

    public StsMovieClass()
    {
        userName = "Movies";
    }

    public void initializeDefaultFields()
    {
        defaultFields = new StsFieldBean[]
            {
            new StsBooleanFieldBean(this, "defaultLoop", "Loop"),
            new StsBooleanFieldBean(this, "defaultCycleVolumes", "Cycle Volumes"),
        };
    }

    public void selected(StsMovie movie)
    {
        super.selected(movie);
        setCurrentObject(movie);
    }

    public StsMovie getCurrentMovie()
    {
        return (StsMovie)currentObject;
    }

    public boolean setCurrentObject(StsObject object)
    {
        return super.setCurrentObject(object);
    }

    public boolean setCurrentMovieName(String name)
    {
        StsMovie newMovie = (StsMovie)getObjectWithName(name);
        return setCurrentObject(newMovie);
    }

    public void close()
    {
        list.forEach("close");
    }
    public boolean getDefaultLoop() { return defaultLoop; }
    public void setDefaultLoop(boolean b)
    {
        if(this.defaultLoop == b) return;
        defaultLoop = b;
//        setDisplayField("defaultLoop", defaultLoop);
    }
    public boolean getDefaultCycleVolumes() { return defaultCycleVolumes; }
    public void setDefaultCycleVolumes(boolean b)
    {
        if(this.defaultCycleVolumes == b) return;
        defaultCycleVolumes = b;
//        setDisplayField("defaultCycleVolumes", defaultCycleVolumes);
    }
}
