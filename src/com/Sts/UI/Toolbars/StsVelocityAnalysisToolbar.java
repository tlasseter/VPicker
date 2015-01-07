package com.Sts.UI.Toolbars;

import com.Sts.Actions.Wizards.VelocityAnalysis.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;

public class StsVelocityAnalysisToolbar extends StsToolbar implements StsSerializable
{
	public static final String NAME = "Velocity Analysis Toolbar";
    public static final boolean defaultFloatable = true;

    /** button filenames (also used as unique identifier button names) */
	public static final String EDIT_VELOCITY = "editVelocity";
//	public static final String VEL_OR_STACK = "displayVelOrStack";
	public static final String DISPLAYSTACK = "displayStack";
	public static final String DISPLAY_TYPE = "displayTypeString";
	public static final String COMPUTE_STACKED_VOLUME = "computeStackedVolume";
	public static final String COMPUTE_SEMBLANCE_VOLUME = "computeSemblanceVolume";
    public static final String UPDATE_VVS = "Update VVS";
    public static final String UNDO = "undo";
	public static final String NMO = "NMO";
	public static final String MODE = "Mode";

	transient StsToggleButton modeToggleButton = null;
    transient StsToggleButton editToggleButton = null;
	transient StsPreStackLineSet currentVolume = null;
	transient StsPreStackVelocityModel velocityModel = null;
	transient StsPreStackLineSetClass lineSetClass = null;
	transient StsToggleButton flattenToggleButton = null;
	transient StsModel model;

	static final long serialVersionUID = 1l;

	public StsVelocityAnalysisToolbar()
	{
        super(NAME);
    }

	public StsVelocityAnalysisToolbar(StsWin3dBase win3d)
	{
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
		this.model = win3d.getModel();
        StsActionManager actionManager = win3d.getActionManager();
        setBorder(BorderFactory.createEtchedBorder());

		editToggleButton = new StsToggleButton(EDIT_VELOCITY, "Add or edit velocity profile", actionManager, StsVelocityAnalysisEdit.class);
		editToggleButton.addIcons("editVelocity", "editVelocityEnd");
		editToggleButton.doClick(); //go ahead and turn picking on - user will want to start picking right away - SWC 5/29/09
		add(editToggleButton);

		flattenToggleButton = new StsToggleButton(NMO, "Select whether you want to flatten or unflatten traces", this, "setFlatten", "setUnFlatten", new Object[0]);
		flattenToggleButton.addIcons("flatten", "unflatten");
		add(flattenToggleButton);

		StsPreStackLineSet currentSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
		setPreStackSeismicSet(currentSet);

		//stackButton = new StsToggleButton(DISPLAYSTACK, "Stack and display current data.", this, DISPLAYSTACK);
		//stackButton.addIcons(DISPLAYSTACK, DISPLAYSTACK);
		//add(stackButton);

		StsComboBoxFieldBean displayTypeCombo = new StsComboBoxFieldBean(lineSetClass, DISPLAY_TYPE, null, StsPreStackLineSetClass.DISPLAY_TYPE_STRINGS);
		displayTypeCombo.setName(DISPLAY_TYPE);
		add(displayTypeCombo);

		StsComboBoxFieldBean stackOptionCombo = new StsComboBoxFieldBean(lineSetClass, "stackOptionString", null, StsPreStackLineSetClass.STACK_STRINGS);
		add(stackOptionCombo);

		StsButton computeStackedVolumeButton = new StsButton(DISPLAYSTACK, "Compute Stacked PostStack3d", lineSetClass, "computeStackedVolume");
		add(computeStackedVolumeButton);

		StsButton computeSemblanceVolumeButton = new StsButton(DISPLAYSTACK, "Compute Semblance PostStack3d", lineSetClass, COMPUTE_SEMBLANCE_VOLUME);
		add(computeSemblanceVolumeButton);
        
        addSeparator();
		addCloseIcon(model.win3d);
		setMinimumSize();
        return true;
    }

	public void setPreStackSeismicSet(StsPreStackLineSet lineSet)
	{
		if (lineSet == null)
		{
			return;
		}
		this.currentVolume = lineSet;
		//lineSetClass = lineSet.lineSetClass;
		lineSetClass = (StsPreStackLineSetClass)StsModel.getCurrentModel().getStsClass(lineSet.getClass());
		velocityModel = lineSet.velocityModel;
		if (flattenToggleButton != null)
			flattenToggleButton.setSelected(lineSetClass.getFlatten());
	}

	public String getProfileOptionString()
	{
		if (lineSetClass == null)
			return null;
		return lineSetClass.getProfileOptionString();
	}

	public void setProfileOptionString(String option)
	{
		if (lineSetClass == null)
			return;
		lineSetClass.setProfileOptionString(option);
	}

    /** If flatten toggle is selected, setFlatten() which returns false if we cannot flatten; if false, programmatically deselect toggle button */
    public void setFlatten()
	{
		if (lineSetClass == null) return;
		if(!lineSetClass.setFlatten())
        {
            flattenToggleButton.setSelected(false);
            new StsMessage(model.win3d, StsMessage.INFO, "No velocity profiles available.\nCannot flatten unNMOed gathers.");
        }
    }
    
    /** If flatten toggle is deselected, setFUnFlatten() which returns false if we cannot unflatten; if false, programmatically select toggle button */
	public void setUnFlatten()
	{
		if (lineSetClass == null) return;
		if(!lineSetClass.setUnFlatten())
        {
           flattenToggleButton.setSelected(true);
           new StsMessage(model.win3d, StsMessage.INFO, "No velocity profiles available.\nCannot unflatten NMOed gathers.");
        }
    }

    /** If a prestackVolume is selected, set flatten toggle button to default position: flatten if volume is NMOed, unflatten otherwise. */
    public void checkSetFlattenToggleButton(boolean isNMOed)
    {
        if(isNMOed)
            lineSetClass.setFlatten();

        else
            lineSetClass.setUnFlatten();

        flattenToggleButton.setSelected(isNMOed);
    }

    /*public void displayStack()
	{
		if (lineSetClass == null)
			return;
		lineSetClass.setDisplayVelocity(!stackButton.isSelected());
	}*/

	static public boolean checkAddToolbar(StsModel model, StsPreStackLineSet lineSet, boolean floatable)
	{
        StsVelocityAnalysisToolbar velocityAnalysisToolbar;

        if (model.win3d == null)
			return false;

		StsPreStackVelocityModel velocityModel = lineSet.velocityModel;
		if (velocityModel == null)
		{
			model.win3d.removeToolbar(NAME);
			return false;
		}
		if(!model.win3d.hasToolbarNamed(NAME))
		{
			velocityAnalysisToolbar = new StsVelocityAnalysisToolbar(model.win3d);
			model.win3d.addToolbar(velocityAnalysisToolbar);
			velocityAnalysisToolbar.setPreStackSeismicSet(lineSet);
		}
        else
            velocityAnalysisToolbar = (StsVelocityAnalysisToolbar)model.win3d.getToolbarNamed(NAME);
        velocityAnalysisToolbar.checkSetFlattenToggleButton(lineSet.getIsNMOed());
        velocityAnalysisToolbar.setVisible(true);
        if ( !velocityAnalysisToolbar.editToggleButton.isSelected() ) velocityAnalysisToolbar.editToggleButton.doClick();
        return true;
    }
}
