package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsSubType implements StsTreeObjectI
{
	StsClass parentClass;
	String name;
	byte type;

	public StsSubType(StsClass parentClass, String name, byte type)
	{
		this.parentClass = parentClass;
		this.name = name;
		this.type = type;
	}
    public void popupPropertyPanel() { return; }     
	public StsFieldBean[] getDisplayFields(){ return null; }
	public StsFieldBean[] getPropertyFields(){ return null; }
	public StsFieldBean[] getDefaultFields(){ return null; }
	public Object[] getChildren(){ return parentClass.getObjectListOfType(type); }
	public StsObjectPanel getObjectPanel(){ return null; }
	public boolean anyDependencies() { return false; }
	public boolean canExport() { return false; }
	public boolean export() { return false; }
	public boolean canLaunch() { return false; }
	public boolean launch() { return false; }
	public String getName() { return name; }
    public void treeObjectSelected() {}
}
