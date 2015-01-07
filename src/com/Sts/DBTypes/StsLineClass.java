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

public class StsLineClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassTimeDepthDisplayable
{
    boolean displayNames = false;
    boolean displayLines = true;

    public StsLineClass()
    {
    }

    public void initializeDisplayFields()
    {
        StsLine.initColors();

        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displayNames", "Names"),
            new StsBooleanFieldBean(this, "displayLines", "Enable")
        };
    }

    public boolean getDisplayNames() { return displayNames; }
    public boolean getDisplayLines() { return displayLines; }

    public void setDisplayNames(boolean display)
    {
        displayNames = display;
        currentModel.win3dDisplayAll();
    }

    public void setDisplayLines(boolean display)
    {
        displayLines = display;
        currentModel.win3dDisplayAll();
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsLine line = (StsLine)iter.next();
            if(displayLines)
                line.display(glPanel3d, displayNames);
        }
    }

    public void displayClass2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsLine line = (StsLine)iter.next();
            if(displayLines)
                line.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }

    public void projectRotationAngleChanged()
    {
        forEach("projectRotationAngleChanged");
    }

	/** classInitialize after db is loaded */
    public boolean dbInitialize()
    {
        initLinesNotOnSections();
        StsSectionClass sectionClass = (StsSectionClass)currentModel.getStsClass(StsSection.class);
        if(sectionClass == null) return true;
        sectionClass.checkSections();
        sectionClass.initSectionsNotOnSections();

        boolean allLinesInitialized = initLinesOnSections();
        boolean allSectionsInitialized = sectionClass.initSections();

        int noIter = 0;
        while (noIter < 3 && (!allLinesInitialized || !allSectionsInitialized))
        {
            allLinesInitialized = initLinesOnSections();
            allSectionsInitialized = sectionClass.initSections();
            noIter++;
        }
        return allLinesInitialized && allSectionsInitialized;
    }

    /** This is the first step in the wells-and-sections initialization process:
     *  classInitialize wells which are not on sections */
    public void initLinesNotOnSections()
    {
        int nLines = getSize();

        for(int n = 0; n < nLines; n++)
        {
            StsLine line = (StsLine)getElement(n);
            if (line != null && line.getOnSection() == null)
                line.initialize();
        }
    }

    /** This is the third step in the wells-and-sections initialization process after
     *  initializing sections not on other sections: classInitialize wells on sections
     *  whether section is classInitialize or not.  We iterate until both wells and sections
     *  are properly initialized.
     */
    public boolean initLinesOnSections()
    {
        boolean allInitialized = true;

        int nLines = getSize();

        for(int n = 0; n < nLines; n++)
        {
            StsLine line = (StsLine)getElement(n);
            if (line != null && line.getOnSection() != null)
                if(!line.initialize()) allInitialized = false;
        }

        return allInitialized;
    }
}
