
/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie Class instantiated by the movie wizard.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsArbitraryLine extends StsMainObject implements StsTreeObjectI
{
    protected StsObjectRefList hinges;

    transient private StsModel model = null;

    static StsObjectPanel objectPanel = null;
    /**
     * Instance of this object
     */
    transient public static StsArbitraryLine aline = null;

    static public final StsFieldBean[] displayFields = null;

    static public final StsFieldBean[] propertyFields = new StsFieldBean[]
    {
        new StsStringFieldBean(StsArbitraryLine.class, "name", true, "Name:"),
        new StsColorComboBoxFieldBean(StsArbitraryLine.class, "stsColor", "Color:")
    };
    /**
     * Default constructor
     */
    public StsArbitraryLine()
    {

    }

    /**
     * constructor
     */
    public StsArbitraryLine(StsModel model)
    {
		super(false);
        this.model = model;
		hinges = StsObjectRefList.constructor(2, 2, "hinges", this);
		addToModel();
//		refreshObjectPanel();
   }

    public boolean initialize(StsModel model)
    {
        try
        {
            this.model = model;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsArbitraryLine.classInitialize() failed.", e, StsException.WARNING);
            return false;
        }
    }

    /**
     * Add a hinge point
     */
    public boolean addHinge(StsLine line)
    {
        hinges.add(line);
        return false;
    }

    /**
     * Arbitrary Line selected on the Object tree
     */
    public void treeObjectSelected()
    {
        getArbitraryLineClass().selected(this);
        currentModel.getGlPanel3d().checkAddView(StsView3d.class);
        currentModel.win3dDisplayAll();
    }

    static public StsArbitraryLineClass getArbitraryLineClass()
    {
        return (StsArbitraryLineClass)currentModel.getCreateStsClass(StsArbitraryLine.class);
    }

    public boolean anyDependencies()
    {
        return false;
    }
    public StsFieldBean[] getDisplayFields() { return displayFields; }
    public StsFieldBean[] getPropertyFields() { return propertyFields; }
    public Object[] getChildren() { return new Object[0]; }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

}
