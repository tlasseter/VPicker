package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;

public class StsVvsUpdateToolbar extends StsToolbar implements StsSerializable
{
	public static final String NAME = "VVS Update Toolbar";
    public static final boolean defaultFloatable = true;

    public static final String UPDATE_VVS = "Update VVS";
    transient StsButton updateVVSButton = null;

	static final long serialVersionUID = 1l;

	public StsVvsUpdateToolbar()
	{
        super(NAME);
    }

	public StsVvsUpdateToolbar(StsWin3dBase win3d)
	{
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        try
        {
            setBorder(BorderFactory.createEtchedBorder());
            updateVVSButton = new StsButton(UPDATE_VVS, "Update any variable velocity stack views in family.", this, "updateVVS", win3d);
            add(updateVVSButton);
            addSeparator();
            if(win3d != null) addCloseIcon(win3d);
            setMinimumSize();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initialize", e);
            return false;
        }
    }

    public void updateVVS(StsWin3dBase win3d)
    {
        win3d.getWindowFamily().viewObjectChangedAndRepaint(this, this);
    }

    public void updateVVS(StsWin3d win3d)
    {
        win3d.getWindowFamily().viewObjectChangedAndRepaint(this, this);
    }

    static public boolean checkAddToolbar(StsWin3dBase win3d)
	{
        if (win3d == null) return false;
        StsVvsUpdateToolbar vvsUpdateToolbar = (StsVvsUpdateToolbar)win3d.getToolbarNamed(NAME);
        if(vvsUpdateToolbar != null)
            vvsUpdateToolbar.setVisible(true);
        else
            win3d.addToolbar(new StsVvsUpdateToolbar(win3d));
        return true;
    }
    static public void main(String[] args)
    {
        try
        {
			javax.swing.UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
            StsVvsUpdateToolbar toolbar = new StsVvsUpdateToolbar(null);
            StsToolkit.createDialog(toolbar);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}