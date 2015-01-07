package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 10, 2007
 * Time: 12:00:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFractureClass extends StsSectionClass implements StsClassDisplayable, StsClassTimeDepthDisplayable
{
    boolean displaySectionEdges = false;

    public StsFractureClass()
    {
        userName = "Fractures";
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        int nFractures = getSize();
        for(int n = 0; n < nFractures; n++)
        {
            StsFracture fracture = (StsFracture)getElement(n);
            fracture.display(glPanel3d, displaySectionEdges);
        }
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
        {
                new StsBooleanFieldBean(this, "displaySections", "Enable")
        };
    }
}
