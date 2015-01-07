package com.Sts.UI.ObjectPanel;

import com.Sts.Interfaces.*;
import com.Sts.Reflect.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: A treeNodeProxy references an StsTreeObjectI which may or may not exist.
 * Until it's created, null or false returns are made from this proxy. </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsTreeNodeProxy implements StsTreeObjectI
{
    Class treeObjectClass;
    StsTreeObjectI[] instanceList = new StsTreeObjectI[1];
    StsTreeObjectI instance = null;
    StsMethod getInstanceMethod;

    public StsTreeNodeProxy(String className)
    {
        try
        {
            treeObjectClass = Class.forName(className);
            getInstanceMethod = new StsMethod(className, "getInstance", new Object[] { instanceList } );
        }
        catch(Exception e)
        {
            StsException.systemError("Failed to construct StsTreeNodeProxy for class " + className);
        }
    }

    private void getInstance()
    {
        getInstanceMethod.invokeStaticMethod(new Object[] { instanceList } );
        instance = instanceList[0];
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(instance == null) getInstance();
        if(instance == null) return null;
        return instance.getDisplayFields();
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(instance == null) getInstance();
        if(instance == null) return null;
        return instance.getPropertyFields();
    }

	public StsFieldBean[] getDefaultFields()
	{
		if(instance == null) getInstance();
		if(instance == null) return null;
		return instance.getDefaultFields();
	}

    public Object[] getChildren()
    {
         if(instance == null) getInstance();
         if(instance == null) return new Object[0];
         return instance.getChildren();
    }

    public StsObjectPanel getObjectPanel()
    {
         if(instance == null) getInstance();
         if(instance == null) return null;
         return instance.getObjectPanel();
    }

    public boolean anyDependencies()
    {
         if(instance == null) getInstance();
         if(instance == null) return false;
         return instance.anyDependencies();
    }
    public boolean export()
    {
         if(instance == null) getInstance();
         if(instance == null) return false;
         return instance.export();
    }
    public boolean canExport()
    {
         if(instance == null) getInstance();
         if(instance == null) return false;
         return instance.canExport();
    }
    public boolean launch()
    {
         if(instance == null) getInstance();
         if(instance == null) return false;
         return instance.launch();
    }
    public boolean canLaunch()
    {
         if(instance == null) getInstance();
         if(instance == null) return false;
         return instance.canLaunch();
    }
    public String getName()
    {
         if(instance == null) getInstance();
         if(instance == null) return null;
         return instance.getName();
    }
    public void treeObjectSelected()
    {
         if(instance == null) getInstance();
         if(instance == null) return;
         instance.treeObjectSelected();
    }
    
    public void popupPropertyPanel()
    {
         if(instance == null) getInstance();
         if(instance == null) return;
         instance.popupPropertyPanel();
    }

}
