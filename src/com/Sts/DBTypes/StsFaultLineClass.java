package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

import java.util.*;

public class StsFaultLineClass extends StsLineClass implements StsSerializable, StsClassDisplayable, StsClassTimeDepthDisplayable
{
    public StsFaultLineClass()
    {
        userName = "Fault Lines";                                                              
    }

    public void initializeDisplayFields()
    {
//        initColors(StsFaultLine.getStaticDisplayFields());

        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displayNames", "Names"),
            new StsBooleanFieldBean(this, "displayLines", "Enable")
        };
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsFaultLine faultLine = (StsFaultLine)iter.next();
            if(displayLines)
                faultLine.display(glPanel3d, displayNames);
        }
    }

    public void displayClass2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsFaultLine faultLine = (StsFaultLine)iter.next();
            if(displayLines)
                faultLine.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }
}
