package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

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
public abstract class StsPanelProperties extends StsSerialize implements StsDialogFace, StsSerializable
{
	String panelTitle;
	String fieldName;
    transient Object parentObject = null;

	transient StsFieldBean[] propertyBeans;
	transient public Object originalProperties;
	private transient boolean changed = false;

	public StsPanelProperties()
	{
	}

	public StsPanelProperties(String panelTitle, String fieldName)
	{
		this.panelTitle = panelTitle;
		this.fieldName = fieldName;
	}

	public StsPanelProperties(Object parentObject, String panelTitle, String fieldName)
	{
		this.parentObject = parentObject;
		this.panelTitle = panelTitle;
		this.fieldName = fieldName;
    }

    abstract void initializeDefaultProperties(Object defaultProperties);

    public void setParentObject(Object parentObject)
    {
        this.parentObject = parentObject;
    }

    public StsDialogFace getEditableCopy()
    {
        StsPanelProperties editableProperties = (StsPanelProperties)StsToolkit.copyAllObjectFields(this);
        editableProperties.parentObject = parentObject;
        editableProperties.setOriginalProperties(this);
        return editableProperties;
    }

    public void setOriginalProperties(StsPanelProperties properties)
    {
        originalProperties = properties;
    }
 /*
    public void saveState()
	{
//		editedProperties = StsToolkit.copyAllObjectFields(this);
	}
*/
	abstract public void initializeBeans();

	public StsFieldBean[] getPropertyBeans()
	{
//		saveState();
		if(propertyBeans == null)
           initializeBeans();
		return propertyBeans;
	}

	/** return true if properties are saved */
	public void dialogSelectionType(int selectionType)
	{
		if(!saveProperties(selectionType)) return;
//		saveState();
		commitChanges();
		update();
		display();
	}

	boolean saveProperties(int selectionType)
	{
		changed = !StsToolkit.compareObjects(this, originalProperties);
		if(!changed) return false;
//        changed = false;
        if(selectionType == StsDialogFace.CLOSING || selectionType == StsDialogFace.CANCEL)
		{
			boolean save = StsYesNoDialog.questionValue(currentModel.win3d, panelTitle+" values have changed; do you wish to save them?");
			if(!save) return false;
        }
        StsToolkit.copyAllObjectNonTransientFields(this, originalProperties);
        return true;
    }

	public void update()
	{
		if(parentObject != null)
		{
			if(parentObject instanceof StsObject) 
				((StsObject)parentObject).objectPropertiesChanged(this);
			currentModel.viewObjectChanged(this, this);
		}
		changed = false;
	}

	public void display()
	{
	    currentModel.viewObjectRepaint(this, this);
	}

	public void commitChanges()
	{
		if (parentObject != null && parentObject instanceof StsObject)
			((StsObject)parentObject).dbFieldChanged(fieldName, this);
	}

	public boolean isChanged()
	{
		return changed;
	}
	
	public Component getPanel()
	{
		return getPanel(false);
	}
	
	public Component getPanel(boolean expandIt)
	{
        return new StsPropertiesPanel(getPropertyBeans(), panelTitle, expandIt);
	}

    public void setChanged(boolean changed)
    {
        this.changed = changed;
    }
}
