package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.MVC.StsProject;
import com.Sts.UI.Beans.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsObjectListPanel;
import com.Sts.UI.StsOkDialog;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

public class StsMainObject extends StsObject
{
    public String name = null;
    protected boolean isVisible = true;
    protected byte type = StsParameters.NONE;

	static final long serialVersionUID = 11;
    public int numberOfElements = 0;

    public StsMainObject()
    {
    }

    public StsMainObject(boolean persistent)
    {
        super(persistent);
    }

	public StsMainObject(boolean persistent, String name)
	{
		super(persistent);
		setName(name);
    }
    public byte getType() { return type; }
    public String getTypeAsString() { return new String("Default-" + type); }
	public void setType(byte type) { this.type = type; }

    public String getUnitsString() { return "none"; }

    public void setName(String name)
    {
        this.name = name;
        if(isPersistent()) dbFieldChanged("name", name);
        if(currentModel != null) currentModel.refreshObjectPanel(this);
    }

    public String getName() { return name; }

    static public String getObjectName(StsMainObject object)
    {
        if(object == null) return "null";
        else return object.getName(); 
    }

    public void setIsVisible(boolean isVisible)
    {
        if(this.isVisible == isVisible) return;
        this.isVisible = isVisible;
        dbFieldChanged("isVisible", isVisible);
        // this is deadly: we might have a whole series of these called on different objects,
        // each of which will trigger a repaint, which is the situation we have tried to eliminate.
        // So in making this call, the caller has the responsibility at the appropriate time
        // to call for a repaint.  In many cases, a number of other things are computed and changed
        // as well before the repaint is called.
        // There are currently 42 calls to this method.  A quick check shows they are doing a lot of
        // stuff before a repaint is called.
        // So I've commented out this line for now.  TJL 2/2/08
        // currentModel.viewObjectRepaint(this);
    }

    public boolean getIsVisible() { return isVisible; }
    public boolean isType(byte type) { return this.type == type; }


    static public Comparator getNameComparator() { return new NameComparator(); }

    static public final class NameComparator implements Comparator
    {
        NameComparator()
        {
        }

        // order by versions and then order alphabetically
        public int compare(Object o1, Object o2)
        {
            StsMainObject so1 = (StsMainObject)o1;
            StsMainObject so2 = (StsMainObject)o2;
            String name1 = so1.getName();
            String name2 = so2.getName();
            return name1.compareTo(name2);
        }
    }

    static public String[] getNamesFromObjects(StsMainObject[] objects)
    {
        if(objects == null) return new String[0];
        int nObjects = objects.length;
        String[] names = new String[nObjects];
        for(int n = 0; n < nObjects; n++)
            names[n] = objects[n].getName();
        return names;
    }

    static public StsMainObject getListObjectWithName(StsMainObject[] objects, String name)
    {
        if(objects == null || name == null) return null;
        int nObjects = objects.length;
        for(int n = 0; n < objects.length; n++)
            if(objects[n].getName().equals(name)) return objects[n];
        return null;
    }
   
    public String toString()
    {
        return getName();
    }

    public String toDebugString()
    {
        return super.toDebugString() + " Name: " + name;
    }

	public StsFieldBean[] getDefaultFields() { return null; }

    public int getNumberOfElements() { return numberOfElements; }
    public void setNumberOfElements(int num) { numberOfElements = num; }

}
