package com.Sts.UI;

import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Oct 3, 2008
 * Time: 8:20:23 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsViewItem implements Serializable
{
    public String name;
    public String shortName;
    public String[] viewClassNames;

    transient public Class[] viewClasses;
    transient public String[] viewNames;
    transient private boolean selected = false;

    public StsViewItem()
    {
    }

    public StsViewItem(Class viewClass)
    {
        name = (String)StsMethod.invokeStaticMethod(viewClass, "getStaticViewName", new Class[0]);
        shortName = (String)StsMethod.invokeStaticMethod(viewClass, "getStaticShortViewName", new Class[0]);
        viewNames = new String[] { name };
        viewClasses = new Class[] { viewClass };
        viewClassNames = new String[] { viewClass.getName() };
    }

    public StsViewItem(String name, Class viewClass)
    {
        this.name = name;
        viewNames = new String[] { name };
        viewClasses = new Class[] { viewClass };
        viewClassNames = new String[] { viewClass.getName() };
    }

    public StsViewItem(String name, Class viewClass, String shortName)
    {
        this.name = name;
        this.shortName = shortName;
        viewNames = new String[] { name };
        viewClasses = new Class[] { viewClass };
        viewClassNames = new String[] { viewClass.getName() };
    }

    public StsViewItem(String name, Class[] viewClasses, String[] viewNames)
    {
        this.name = name;
        this.shortName = name;
        this.viewClasses = viewClasses;
        this.viewNames = viewNames;
        int nViewClasses = viewClasses.length;
        viewClassNames = new String[nViewClasses];
        for(int n = 0; n < nViewClasses; n++)
            viewClassNames[n] = viewClasses[n].getName();
    }

    public String toString() { return name; }

    public Class[] getViewClasses() { return viewClasses; }

    public Class getViewClass() { return viewClasses[0]; }
    /*
    public void selected()
    {
        try
        {
            selectedMethod.invokeInstanceMethod(win3d, new Object[]{viewClasses, name});
        }
        catch(Exception e)
        {
            StsException.systemError(this, "selected");
        }
    }
    */
    public String getName() { return name; }

    /** assume each if these viewItems is for a single view;
     *  return the viewClasses list for these views.
     * @param viewItems
     * @return viewClasses
     */
    static public Class[] getViewClasses(Object[] viewItems)
    {
        int nViewItems = viewItems.length;
        Class[] viewClasses = new Class[nViewItems];
        for(int n = 0; n < nViewItems; n++)
        {
            StsViewItem viewItem = (StsViewItem)viewItems[n];
            viewClasses[n] = viewItem.getViewClass();
        }
        return viewClasses;
    }

    /** return list of viewNames for these view items.
     *
     * @param viewItems
     * @return list of viewItem names
     */
    static public String[] getViewNames(Object[] viewItems)
    {
        int nViewItems = viewItems.length;
        String[] viewNames = new String[nViewItems];
        for(int n = 0; n < nViewItems; n++)
            viewNames[n] = ((StsViewItem)viewItems[n]).name;
        return viewNames;
    }

    /** Given these single view items and a set of viewNames for a multiViewPanel,
     *  return the Classes for the desired multiViewPanel
     * @param viewItems
     * @param viewNames
     * @return viewClasses
     */
    static public Class[] getViewClasses(Object[] viewItems, String[] viewNames)
    {
        int nViewNames = viewNames.length;
        Class[] viewClasses = new Class[nViewNames];
        int nFound = 0;
        int nViewItems = viewItems.length;
        for(int n = 0; n < nViewNames; n++)
        {
            String viewName = viewNames[n];
            for(int i = 0; i < nViewItems; i++)
            {
                StsViewItem viewItem = (StsViewItem)viewItems[i];
                if(viewItem.name.equals(viewName))
                {
                    viewClasses[nFound++] = viewItem.getViewClass();
                }
            }
        }
        return viewClasses;
    }

    public boolean getSelected()
    {
        return selected;
    }

    public void setSelected(boolean selected)
    {
        this.selected = selected;
    }

    public void initializeViewClasses()
    {
        int nViewClasses = viewClassNames.length;
        viewClasses = new Class[nViewClasses];
        viewNames = new String[nViewClasses];
        for(int n = 0; n < nViewClasses; n++)
        {
            try
            {
                viewClasses[n] = Class.forName(viewClassNames[n]);
                viewNames[n] = (String)StsMethod.invokeStaticMethod(viewClasses[n], "getStaticViewName", new Class[0]);
            }
            catch(Exception e)
            {
                StsException.outputWarningException(this, "initializeViewClasses", e);
                return;
            }
        }
    }
}

