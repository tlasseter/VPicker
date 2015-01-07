package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>e
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsPlatformOrWellPanel extends StsFieldBeanPanel
{
	boolean isNewPlanSet = true;
	StsWellPlanSet[] wellPlanSets = new StsWellPlanSet[0];
	StsPlatform[] platforms = new StsPlatform[0];
	StsWell nullWell =  new StsWell("none", false);
	StsWell[] availableWells = new StsWell[] { nullWell };
	StsWellPlanSet selectedPlanSet = null;
	String wellType = PLATFORM_WELL;
	String wellPlanSetName = "WellPlan-" + planSetNumber++;
	boolean isNewPlatform = true;
	String platformName = "Platform-" + platformNumber++;
	StsWell drillingWell = null;
	StsPlatform existingPlatform = null;

    private StsDefineWellPlanWizard wizard;

    StsGroupBox newOrOldPlanBox;
	StsBooleanFieldBean newPlanCheckboxBean = new StsBooleanFieldBean(this, "isNewPlan", true, "Define new set of well plans");
	// one of these two will be isVisible depending on whether plan is new or old
    StsStringFieldBean planNameBean = new StsStringFieldBean(this, "wellPlanSetName", "Plan name");
	StsComboBoxFieldBean selectPlanSetBean = new StsComboBoxFieldBean(this, "selectedPlanSet", "Select Plan");

	StsComboBoxFieldBean wellTypeBean = new StsComboBoxFieldBean(this, "wellType", "Specify well type", wellTypeStrings);

	StsBooleanFieldBean newPlatformCheckboxBean = new StsBooleanFieldBean(this, "isNewPlatform", true, "New platform");
	// one of these two will be isVisible depending on whether new or existing platform specified above
	StsStringFieldBean platformNameBean = new StsStringFieldBean(this, "platformName", "Platform name");
	StsComboBoxFieldBean existingPlatformsBean = new StsComboBoxFieldBean(this, "existingPlatform", "Existing platform");

	StsGroupBox drillingWellBox;
	StsComboBoxFieldBean drillingWellsBean = new StsComboBoxFieldBean(this, "drillingWell", "Select Drilling Well");

	public static String PLATFORM_WELL = "Platform well";
    public static String SINGLE_WELL = "Single well";
    public static String PATTERN_WELL = "Pattern well";

	public static int planSetNumber = 0;
	public static int platformNumber = 0;

    public static String[] wellTypeStrings = new String[] { PLATFORM_WELL, SINGLE_WELL, PATTERN_WELL };

//    static String DRILLING_WELL = "Extend drilling well plan";
//    static String[] extendWellPlanOptions = new String[] { DRILLING_WELL };

    public StsPlatformOrWellPanel(StsWizard wizard)
    {
        this.wizard = (StsDefineWellPlanWizard)wizard;
        try
        {
//            layoutPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
		StsModel model = wizard.getModel();
		wellPlanSets = (StsWellPlanSet[])model.getCastObjectList(StsWellPlanSet.class);
		selectPlanSetBean.setListItems(wellPlanSets);
		selectedPlanSet = (StsWellPlanSet)model.getCurrentObject(StsWellPlanSet.class);
		platforms = (StsPlatform[])model.getCastObjectList(StsPlatform.class);
		existingPlatformsBean.setListItems(platforms);
		existingPlatform = platforms.length == 0 ? null: platforms[0];
		StsWell[] wells = (StsWell[])model.getCastObjectList(StsWell.class);
        wells = (StsWell[]) StsMath.arrayAddArray(wells, model.getCastObjectList(StsLiveWell.class));
		availableWells = (StsWell[])StsMath.arrayAddArray(availableWells, wells);
		drillingWellsBean.setListItems(availableWells);
		if(wellPlanSets.length < 1) isNewPlanSet = true;
	    layoutPanel();
    }

    void layoutPanel()
    {
		gbc.fill = GridBagConstraints.NONE;
	    newOrOldPlanBox = new StsGroupBox("1. Select/define well plan.");
		drillingWellBox = new StsGroupBox("2. Select drilling will associated with this plan (already loaded).");

		if(wellPlanSets.length > 0)
			newOrOldPlanBox.add(newPlanCheckboxBean);
		if(isNewPlanSet)
		{
			newOrOldPlanBox.add(planNameBean);
			newOrOldPlanBox.add(wellTypeBean);
			if(wellType == PLATFORM_WELL && platforms.length > 0)
				newOrOldPlanBox.add(newPlatformCheckboxBean);
		}
		else
		{
			newOrOldPlanBox.add(selectPlanSetBean);
		}
		add(newOrOldPlanBox);

		if(wellType == PLATFORM_WELL)
		{
			if (isNewPlatform)
				newOrOldPlanBox.add(platformNameBean);
			else
				newOrOldPlanBox.add(existingPlatformsBean);
		}
	    drillingWellBox.add(drillingWellsBean);
		add(drillingWellBox);
    }

    public String getWellType() { return wellType; }
    public void setWellType(String wellType)
    {
		if (this.wellType == wellType)return;
		// if choice is changed from or to PLATFORM_WELL, we need to layout panel
		boolean rebuild = this.wellType == PLATFORM_WELL || wellType == PLATFORM_WELL;
		this.wellType = wellType;
		if(rebuild) rebuild();
    }

	private void rebuild()
	{
		removeAll();
		repaint();
		layoutPanel();
		wizard.rebuild();
	}

	public boolean getIsNewPlan() { return isNewPlanSet; }
	public void setIsNewPlan(boolean isNewPlan)
	{
		if(this.isNewPlanSet == isNewPlan) return;
		this.isNewPlanSet = isNewPlan;
		if(isNewPlan) selectedPlanSet = null;
		rebuild();
	}

	public void setSelectedPlanSet(StsWellPlanSet planSet)
	{
		selectedPlanSet = planSet;
		planNameBean.setValue(selectedPlanSet.getName());
		StsPlatformClass pc = (StsPlatformClass)wizard.getModel().getStsClass(StsPlatform.class);
		existingPlatform = pc.getWellPlatform(wellPlanSetName);

	}
	public StsWellPlanSet getSelectedPlanSet() { return selectedPlanSet; }

	public boolean getIsNewPlatform() { return isNewPlatform; }
	public void setIsNewPlatform(boolean isNewPlatform)
	{
		if(this.isNewPlatform == isNewPlatform) return;
		this.isNewPlatform = isNewPlatform;
		if(isNewPlatform) resetPlatformName();
		rebuild();
	}

	private void resetPlatformName()
	{
		platformName = "Platform-" + platformNumber++;
	}

	private void resetWellPlanSetName()
	{
		wellPlanSetName = "WellPlan-" + planSetNumber++;
	}

	public void setDrillingWell(StsWell well)
	{
		if(well == nullWell)
			drillingWell = null;
		else
			this.drillingWell = well;
	}
	public StsWell getDrillingWell() { return drillingWell; }

	public void setExistingPlatform(StsPlatform platform)
	{
		this.existingPlatform = platform;
		platformNameBean.setValue(existingPlatform.getName());
	}
	public StsPlatform getExistingPlatform() { return existingPlatform; }

	public String getPlatformName() { return platformName; }
	public void setPlatformName(String name) { this.platformName = name; }

	public String getWellPlanSetName()
	{
		if(this.selectedPlanSet != null)
			return selectedPlanSet.getName();
		else
			return "WellPlanSet-" + planSetNumber;
	}

	public void setWellPlanSetName(String name) { this.wellPlanSetName = name; }

	public StsWellPlanSet getWellPlanSet()
	{
		if(isNewPlanSet)
		{
			selectedPlanSet = new StsWellPlanSet(wellPlanSetName);
			planSetNumber++;
		}
		StsWellPlan currentWellPlan = selectedPlanSet.getCurrentWellPlan();
		if(currentWellPlan == null)
		{
			StsException.systemError("StsPloatformOrWellPanel.getWellPlanSet() failed. WellPlanSet has no well plans.");
			return null;
		}
	/*
		if(drillingWell != null)
		{
			currentWellPlan.setDrillingWell(drillingWell);
		}
	*/
//		selectedPlanSet.addWellPlan(currentWellPlan);
		return selectedPlanSet;
	}

	public StsPlatform getCreatePlatform()
	{
		if(wellType != PLATFORM_WELL) return null;
		if(isNewPlatform)
		{
			StsPlatform platform = new StsPlatform(true);
			platform.setName(platformName);
			return platform;
		}
		else
		{
			return existingPlatform;
		}
    }
}
