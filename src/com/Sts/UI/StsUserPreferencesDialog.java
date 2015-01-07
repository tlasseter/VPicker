package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */


public class StsUserPreferencesDialog extends JDialog implements StsDialogFace
{
    StsBooleanFieldBean variableBean = null;
    StsBooleanFieldBean windowsBean = null;
    boolean saveVariables = true, saveWindows = false;
    StsModel model;
    StsFieldBean[] allBeans = null;

	public StsUserPreferencesDialog(StsModel model, boolean modal)
	{
        this.model = model;
        initializeBeans();
        StsOkCancelDialog dialog = new StsOkCancelDialog(model.win3d, this, "User Preferences", modal);
	}

    private void initializeBeans()
    {
        StsFieldBean[] displayBeans, defaultBeans;
        String name;

        StsProject project = model.getProject();
        displayBeans = project.getDisplayFields();
        defaultBeans = project.getDefaultFields();
        addBeans(displayBeans, defaultBeans);

        StsClass[] classes = model.classList;
        for(int n = 0; n < classes.length; n++)
        {
            StsClass stsClass = classes[n];
            displayBeans = stsClass.getDisplayFields();
            defaultBeans = stsClass.getDefaultFields();
            addBeans(displayBeans, defaultBeans);
        }
    }

    private void addBeans(StsFieldBean[] displayBeans, StsFieldBean[] defaultBeans)
    {
        if(displayBeans == null && defaultBeans == null)
            return;

        if (displayBeans != null)
            allBeans = (StsFieldBean[])StsMath.arrayAddArray(allBeans, displayBeans);

        if (defaultBeans != null)
            allBeans = (StsFieldBean[])StsMath.arrayAddArray(allBeans, defaultBeans);
    }
    
    public Component getPanel(boolean val) { return getPanel(); }
    public Component getPanel()
    {
		StsJPanel dialogPanel = StsJPanel.addInsets();
        StsGroupBox beanPanel = new StsGroupBox("User Preferences");
        beanPanel.gbc.anchor = beanPanel.gbc.WEST;
        beanPanel.gbc.fill = beanPanel.gbc.BOTH;
        dialogPanel.add(beanPanel);

        variableBean = new StsBooleanFieldBean(this, "saveVariables", "Save Variable Settings:");
        windowsBean = new StsBooleanFieldBean(this, "saveWindows", "Save Window Configuration:");
        beanPanel.add(variableBean);
        beanPanel.add(windowsBean);
		return dialogPanel;
    }

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }
    
    public void dialogSelectionType(int type)
    {
        Properties defaultProperties = new Properties();
		try
		{
            if(type == StsDialogFace.OK)
            {
                if(saveVariables)
                {
                    for (int n = 0; n < allBeans.length; n++)
                    {
                        String beanKey = allBeans[n].getBeanKey();
                        if (beanKey == null)continue;
                        Object valueObject = allBeans[n].getValueObject();
                        if (valueObject == null)continue;
                        String objectString = allBeans[n].toString();
                        defaultProperties.put(beanKey, objectString);
                    }
                    String directory = System.getProperty("user.home");
                    StsFile file = StsFile.constructor(directory, "user.preferences");
                    defaultProperties.store(file.getOutputStream(), "Default properties");
                }
                if(saveWindows)
                    model.outputWindowGeometry();
            }
		}
		catch(Exception e)
		{
			StsException.outputException("StsUserPreferencesDialog.dialogSelectionType(OK) failed.", e, StsException.WARNING);
		}
	}
    public void setSaveWindows(boolean val) { saveWindows = true; }
    public void setSaveVariables(boolean val) { saveVariables = true; }
    public boolean getSaveWindows() { return saveWindows; }
    public boolean getSaveVariables() { return saveVariables; }
}
