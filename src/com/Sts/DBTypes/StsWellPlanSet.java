package com.Sts.DBTypes;

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

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsWellPlanSet extends StsMainTimeObject implements StsTreeObjectI, StsInstance3dDisplayable
{
	protected StsObjectRefList wellPlans;
	protected StsWellPlan currentPlan = null;

	public boolean showAll = false;
	public boolean drawVertices = false;
	public boolean displayNames = false;
	public boolean displayPrevious = false;
	public boolean drawCurtain = false;
	public String drawLabelString = StsWell.NO_LABEL;
    public float labelInterval = 100.0f;

    transient double xOrigin;
	transient double yOrigin;

	static protected StsObjectPanel objectPanel = null;

	static final long serialVersionUID = 1L;

	static StsComboBoxFieldBean plansListBean;
    static public StsFieldBean[] displayFields = null;

	public StsWellPlanSet()
	{
	}

	public StsWellPlanSet(String name)
	{
		setName(name);
		initializePlans();
	}

	public boolean initialize(StsModel model)
	{
        resetWellPlanColors();
        if(drawCurtain)
        {
            if(currentPlan != null)
                currentPlan.createSeismicCurtain();
        }
        return true;
	}

	private void initializePlans()
	{
		wellPlans = StsObjectRefList.constructor(2, 2, "wellPlans", this);
		currentPlan = new StsWellPlan(false);
//		wellPlans.add(currentPlan);
	}

	public StsWellPlan getCurrentWellPlan() { return currentPlan; }

	public boolean delete()
	{
		if (!super.delete())
		{
			return false;
		}
		if (currentPlan != null) currentPlan.delete();
		StsObjectRefList.deleteAll(wellPlans);
		checkDeleteSeismicCurtain();
        deleteFromPlatform();
        currentPlan = null;
		return true;
	}

    public void deleteFromPlatform()
    {
        StsPlatformClass platformClass = (StsPlatformClass)currentModel.getStsClass(StsPlatform.class);
        platformClass.deleteWellFromPlatform(getName());
    }

    public StsWellPlan copyLastWellPlan()
	{
		try
		{
			if (wellPlans == null)
			{
				return null;
			}
			StsWellPlan oldPlan = (StsWellPlan) wellPlans.getLast();

			StsWellPlan newPlan = (StsWellPlan)StsToolkit.deepCopy(oldPlan);
			newPlan.setPrevWellPlan(oldPlan);
			newPlan.setIndex(-1); // we need to add this wellPlan to the model later when finished
			newPlan.addToModel();
            newPlan.setStsColor(StsColor.BLUE);
/*
			StsPoint[] oldPoints = oldPlan.getPoints();
			int nPoints = oldPoints.length;
			StsPoint[] newPoints = new StsPoint[nPoints];
			for (int n = 0; n < nPoints; n++)
			{
				newPoints[n] = new StsPoint(oldPoints[n]);
			}
			newPlan.setPoints(newPoints);
*/
			addWellPlan(newPlan);
			return newPlan;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.copyLastPlannedWell() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public StsWellPlan getLastWellPlan()
	{
		if(wellPlans == null) return null;
		return (StsWellPlan)wellPlans.getLast();
	}

	public void addWellPlan(StsWellPlan wellPlan)
	{
		wellPlans.add(wellPlan);
		wellPlan.addTimeStampName();
        currentModel.getProject().setProjectTimeToCurrentTime();  // Forces time to current time so plan displays
		wellPlan.setZDomainSupported(StsLine.TD_TIME_DEPTH);
		setCurrentPlan(wellPlan);
        resetWellPlanColors();
        StsObjectPanel objectPanel = getObjectPanel();
		if (objectPanel.getPanelObject() == this)
		{
            plansListBean.addItem(wellPlan);
            plansListBean.setSelectedItem(wellPlan);
            objectPanel.refreshProperties();
		}
	}

    /** A hack for now to set the colors. These colors aren't persisted so this method is called
     *  when this planSet is loaded from the database.
     */
    private void resetWellPlanColors()
    {
        int nPlans = wellPlans.getSize();
        if(nPlans == 0) return;
        StsWellPlan originalPlan = (StsWellPlan)wellPlans.getElement(0);
        originalPlan.stsColor = new StsColor(StsColor.RED);
        if(nPlans == 1) return;
        StsWellPlan lastPlan = (StsWellPlan)wellPlans.getLast();
        lastPlan.stsColor = new StsColor(StsColor.BLUE);
        if(nPlans == 2) return;
        StsWellPlan prevPlan = (StsWellPlan)wellPlans.getSecondToLast();
        prevPlan.stsColor = new StsColor(StsColor.CYAN);
        if(nPlans == 3) return;
        for(int n = 1; n < nPlans-2; n++)
        {
            StsWellPlan betweenPlan = (StsWellPlan)wellPlans.getElement(n);
            betweenPlan.stsColor = new StsColor(StsColor.GREY);
        }
    }

    public Object[] getWellPlans()
	{
		return wellPlans.getTrimmedList();
    }
	private void setPlansListBeanList()
	{
        StsWellPlan saveCurrentPlan = currentPlan;
        plansListBean.removeAll();
        if(wellPlans == null) return;
        int nWellPlans = wellPlans.getSize();
        for (int n = 0; n < nWellPlans; n++)
		{
			plansListBean.addItem(wellPlans.getElement(n));
		}
        currentPlan = saveCurrentPlan;
        plansListBean.setValueObject(currentPlan);
    }

	/** returns true if changed */
	public boolean setIsVisibleNoDisplay(boolean isVisible)
	{
		if (this.isVisible == isVisible)
		{
			return false;
		}
		this.isVisible = isVisible;
		return true;
    }


	public String getDrawLabelString()
	{
		return drawLabelString;
	}

	public void setDrawLabelString(String labelString)
	{
		drawLabelString = labelString;
        dbFieldChanged("drawLabelString", drawLabelString);
		currentModel.win3dDisplay();
		return;
    }

	public float getLabelInterval()
	{
		return labelInterval;
	}

	public void setLabelInterval(float value)
	{
		labelInterval = value;
        dbFieldChanged("labelInterval", labelInterval);
		currentModel.win3dDisplay();
		return;
    }

	public double getXOrigin()
	{
		if(currentPlan != null)
			return currentPlan.getXOrigin();
		else
			return 0.0;
	}

	public double getYOrigin()
	{
		if(currentPlan != null)
			return currentPlan.getYOrigin();
		else
			return 0.0;
	}

    public void setIsVisible(boolean isVisible)
    {
        super.setIsVisible(isVisible);
        currentModel.win3dDisplayAll();
    }

	public void display(StsGLPanel3d glPanel3d)
	{
 //       long projectTime = currentModel.getProject().getProjectTime();
 //       boolean priorToCurrentTime = true;
		if (showAll)
		{
			for (int i = 0; i < wellPlans.getSize(); i++)
			{
                StsWellPlan plan = (StsWellPlan) wellPlans.getElement(i);
//                if(StsWellPlanSetClass.getWellPlanSetClass().getEnableTime())
//                {
//                   if(plan.getTimeStamp() <= projectTime)
//                      plan.display(glPanel3d, displayNames, drawLabelString, labelInterval);
//                   else
//                      priorToCurrentTime = false;
//                }
//                else
                    plan.display(glPanel3d, displayNames, drawLabelString, labelInterval);
            }
		}
		else
		{
            // Set the current plan to the one closest to the project time
//            if(StsWellPlanSetClass.getWellPlanSetClass().getEnableTime())
//            {
//                if(wellPlans.getSize() > 0)
//                	priorToCurrentTime = false;
//                for (int i = 0; i < wellPlans.getSize(); i++)
//                {
//                    StsWellPlan plan = (StsWellPlan) wellPlans.getElement(i);
//                    if(plan.getTimeStamp() <= projectTime)
//                    {
//                        currentPlan = plan;
//                        priorToCurrentTime = true;
//                    }
//                    else
//                        break;
//                }
//            }
//			if((currentPlan != null) && (priorToCurrentTime))
			if(currentPlan != null)
			{				
                currentPlan.display(glPanel3d, displayNames, drawLabelString, labelInterval);
			}
//			if((displayPrevious) && (priorToCurrentTime))
			if(displayPrevious)
			{				
				int currentIdx = wellPlans.getIndex(currentPlan);
				if ((currentIdx - 1) >= 0)
				{
					((StsWellPlan) wellPlans.getElement(currentIdx - 1)).display(glPanel3d, false, drawLabelString, labelInterval);
				}
			}
		}
//		if((drawCurtain) && (priorToCurrentTime))
		if(drawCurtain)
		{
			if(this.currentPlan == null) return;
			currentPlan.displaySeismicCurtain(glPanel3d);
		}
    }
	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields == null)
        {
            plansListBean = new StsComboBoxFieldBean(StsWellPlanSet.class, "plan", "Well Plans");
            displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsWellPlanSet.class, "name", false, "Plan name"),
                plansListBean,
                new StsBooleanFieldBean(StsWellPlanSet.class, "isVisible", "Enable"),
                new StsBooleanFieldBean(StsWellPlanSet.class, "showAll", "Show All"),
                new StsColorComboBoxFieldBean(StsWellPlanSet.class, "stsColor", "Color:", currentModel.getSpectrum("Basic").getStsColors()),
                new StsBooleanFieldBean(StsWellPlanSet.class, "drawVertices", "Show Vertices"),
                new StsBooleanFieldBean(StsWellPlanSet.class, "displayNames", "Show Names"),
                new StsComboBoxFieldBean(StsWellPlanSet.class, "drawLabelString", "Label Type:", StsWell.LABEL_STRINGS),
                new StsFloatFieldBean(StsWellPlanSet.class, "labelInterval", true, "Label Interval:"),
                new StsBooleanFieldBean(StsWellPlanSet.class, "drawCurtain", "Curtain:"),
        //		new StsColorListFieldBean(StsWellPlanSet.class, "stsColor", "Color"),
                new StsDoubleFieldBean(StsWellPlanSet.class, "xOrigin", false, "X Origin"),
                new StsDoubleFieldBean(StsWellPlanSet.class, "yOrigin", false, "Y Origin")
            };
        }
        setPlansListBeanList();
		return displayFields;
    }

	public void setDrawCurtain(boolean curtain)
	{
		if (drawCurtain == curtain)
		{
			return;
		}
		int nWellPlans = wellPlans.getSize();
	/*
		if (nWellPlans <= 0)
		{
			drawCurtainBean.setValue(drawCurtain);
		}
	*/
		drawCurtain = curtain;
        dbFieldChanged("drawCurtain", drawCurtain);
		this.currentPlan.setDrawCurtain(drawCurtain);
	}
	public boolean getDrawCurtain() { return drawCurtain; }

	public void setDrawVertices(boolean drawVertices) { this.drawVertices = drawVertices; }
	public boolean getDrawVertices() { return drawVertices; }

	public void display2d(StsGLPanel3d glPanel3d, int dir, float coor, boolean axisFlip, boolean xRev, boolean yRev)
	{
//        long projectTime = currentModel.getProject().getProjectTime();
//        boolean priorToCurrentTime = true;
		if (showAll)
		{
			for (int i = 0; i < wellPlans.getSize(); i++)
			{
               StsWellPlan plan = (StsWellPlan) wellPlans.getElement(i);
//               if(StsWellPlanSetClass.getWellPlanSetClass().getEnableTime())
//               {
//                   if(plan.getTimeStamp() <= projectTime)
//                       plan.display2d(glPanel3d, displayNames, dirNo, coor, axisFlip, xRev, yRev);
//                   else
//                       priorToCurrentTime = false;
//               }
//               else
                   plan.display2d(glPanel3d, displayNames, dir, coor, axisFlip, xRev, yRev);
			}
		}
		else
		{
            // Set the current plan to the one closest to the project time
//            if(StsWellPlanSetClass.getWellPlanSetClass().getEnableTime())
//            {
//                priorToCurrentTime = false;
//                for (int i = 0; i < wellPlans.getSize(); i++)
//                {
//                    StsWellPlan plan = (StsWellPlan) wellPlans.getElement(i);
//                    if(plan.getTimeStamp() <= projectTime)
//                    {
//                        currentPlan = plan;
//                        priorToCurrentTime = true;
//                    }
//                    else
//                        break;
//                }
//            }
//			if((currentPlan != null) && (priorToCurrentTime))			
			if(currentPlan != null)
			{
				currentPlan.display2d(glPanel3d, displayNames, dir, coor, axisFlip, xRev, yRev);
			}
//			if((displayPrevious) && (priorToCurrentTime))			
			if(displayPrevious)
			{
				int currentIdx = wellPlans.getIndex(currentPlan);
				if ( (currentIdx - 1) >= 0)
				{
					( (StsWellPlan) wellPlans.getElement(currentIdx - 1)).display2d(glPanel3d, displayNames, dir, coor, axisFlip, xRev, yRev);
				}
			}
		}
    }

	public Object getPlan()
	{
		return currentPlan;
	}

	public void setPlan(Object plannedWell)
	{
		if(plannedWell == currentPlan) return;
		setCurrentPlan((StsWellPlan) plannedWell);
//		setBeachballColors(currentPlan.getColor());
//		setShowVertices(currentPlan.getDrawVertices());
		currentModel.win3dDisplayAll();
    }

	public void setCurrentPlan(StsWellPlan wellPlan)
	{
	    checkDeleteSeismicCurtain();
		this.currentPlan = wellPlan;
		if(drawCurtain) wellPlan.createSeismicCurtain();
	}

	private void checkDeleteSeismicCurtain()
	{
		if(currentPlan != null)
        {
            currentPlan.deleteSeismicCurtain();
            currentPlan.setDrawCurtain(false);
        }
    }

	public boolean getDisplayNames()
	{
		return displayNames;
	}

	public void setDisplayNames(boolean displayNames)
	{
		this.displayNames = displayNames;
        dbFieldChanged("displayNames", displayNames);
		currentModel.win3dDisplayAll();
	}

	public boolean getShowAll()
	{
		return showAll;
	}

	public void setShowAll(boolean showAll)
	{
		this.showAll = showAll;
        dbFieldChanged("showAll", showAll);
		currentModel.win3dDisplayAll();
	}

	public boolean getDisplayPrevious()
	{
		return displayPrevious;
	}

	public void setDisplayPrevious(boolean previous)
	{
		this.displayPrevious = previous;
        dbFieldChanged("displayPrevious", displayPrevious);
		currentModel.win3dDisplayAll();
    }
	public void setPreviousColor(StsColor color)
	{
		int currentIdx = wellPlans.getIndex(currentPlan);
		if ( (currentIdx - 1) >= 0)
		{
			wellPlans.getElement(currentIdx - 1).setStsColor(color);
		}
    }

    public StsColor getStsColor()
    {
        if(this.currentPlan == null) return null;
        return currentPlan.getStsColor();
    }

    public void setStsColor(StsColor color)
    {
        if(this.currentPlan == null) return;
        currentPlan.setStsColor(color);
    }

    public StsFieldBean[] getPropertyFields()
	{
		return null;
	}

	public Object[] getChildren()
	{
		return new Object[0];
	}

	public boolean anyDependencies()
	{
		return false;
    }

	public void treeObjectSelected()
	{
		currentModel.getStsClass(StsWellPlanSet.class).selected(this);
	}
/*
	public boolean export()
	{
		if (this.currentPlan == null)
		{
			return false;
		}
		return currentPlan.export(StsParameters.TD_TIME_DEPTH_STRING);
    }
*/
    public boolean canExport() { return true; }
	public boolean export()
	{
		if (this.currentPlan == null) return false;
		boolean ok = StsYesNoDialog.questionValue(currentModel.win3d, "Export complete well plan? Otherwise only the path will be exported.");
	    if(ok)
			return currentPlan.exportPlan(this.getName());
		else
			return currentPlan.exportPlanWell(StsParameters.TD_TIME_DEPTH_STRING);
	}

	public StsObjectPanel getObjectPanel()
	 {
		 if (objectPanel == null)
		 {
			 objectPanel = StsObjectPanel.constructor(this, true);
		 }
		 return objectPanel;
    }
}
