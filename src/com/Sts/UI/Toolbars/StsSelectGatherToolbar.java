package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;


public class StsSelectGatherToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Select Gather Toolbar";
    public static final boolean defaultFloatable = true;

    public static final String PREVIOUS = "moviePrevious";
    public static final String NEXT = "movieNext";

    transient StsWin3dBase window;

	static final long serialVersionUID = 1l;

	public StsSelectGatherToolbar()
	{
        super(NAME);
    }

	public StsSelectGatherToolbar(StsWin3dBase win3d)
	{
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.window = win3d;
		setBorder(BorderFactory.createEtchedBorder());

        add(new StsButton(PREVIOUS, "Go to previously picked velocity function (PG_DOWN)", this, "previousProfile"));
        add(new StsButton(NEXT, "Go to next velocity function to pick (PG_UP)", this, "nextProfile"));
        StsComboBoxFieldBean profileOptionCombo = new StsComboBoxFieldBean(this, "profileOptionString", null, StsPreStackLineSet3dClass.PROFILE_STRINGS);
        add(profileOptionCombo);

        addSeparator();
		//addCloseIcon(model.win3d);
		setMinimumSize();
        return true;
    }

    // All these methods need to query for the current prestack seismic object since the user may have changed it. Therefore
    // cannot use global currentVolume and lineSetClass.
    public String getProfileOptionString()
    {
        StsPreStackLineSet currentVolume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if(currentVolume == null) return null;
        return currentVolume.lineSetClass.getProfileOptionString();
    }
	public void setProfileOptionString(String option)
	{
        StsPreStackLineSet currentVolume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if(currentVolume == null) return;
        currentVolume.lineSetClass.setProfileOptionString(option);
    }

    public void previousProfile()
    {
		StsPreStackLineSet currentVolume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
		if(currentVolume != null) currentVolume.previousProfile(window);
//        window.resetFamilyViews();
    }

    public void nextProfile()
    {
        StsPreStackLineSet currentVolume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (currentVolume != null) currentVolume.nextProfile(window);
//        window.resetFamilyViews();
    }

}
