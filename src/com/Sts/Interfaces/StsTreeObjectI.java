package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;

public interface  StsTreeObjectI
{
    public StsFieldBean[] getDisplayFields();
    public StsFieldBean[] getPropertyFields();
	public StsFieldBean[] getDefaultFields();
	public Object[] getChildren();
    public StsObjectPanel getObjectPanel();
    public boolean anyDependencies();
    public boolean canExport();
    public boolean export();
    public boolean canLaunch();
    public boolean launch();
    public String getName();
    public void treeObjectSelected();
    public void popupPropertyPanel();
}
