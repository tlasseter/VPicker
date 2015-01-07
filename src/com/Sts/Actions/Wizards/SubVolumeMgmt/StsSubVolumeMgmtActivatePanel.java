package com.Sts.Actions.Wizards.SubVolumeMgmt;

import com.Sts.Actions.Wizards.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSubVolumeMgmtActivatePanel extends StsJPanel
{
    private StsSubVolumeMgmtWizard wizard;
    private StsSubVolumeMgmtActivate wizardStep;

    transient public StsObjectTree objectTreePanel;
    StsTreeNode subVolumeNode = null;

    public StsSubVolumeMgmtActivatePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsSubVolumeMgmtWizard)wizard;
        this.wizardStep = (StsSubVolumeMgmtActivate)wizardStep;

        try
        {
            constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void constructBeans()
    {
        StsModel model = wizard.getModel();
        Dimension panelSize = new Dimension(300, 400);
        objectTreePanel = new StsObjectTree(model, panelSize);
        objectTreePanel.setRootNode(model.getWorkflowPlugIn().getSubVolumeNode());
        objectTreePanel.finalizeTreeModel();
        objectTreePanel.refreshTree();
    }

    public void initialize()
    {

    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.BOTH;
        add(objectTreePanel);
    }

    // Methods for StsTreeObjectI interface
    public boolean canExport() { return false; }
    public boolean export() { return false; }
    public boolean canLaunch() { return false; }
    public boolean launch() { return false; }
    public String getName() { return null; }
    public StsFieldBean[] getDisplayFields() { return null; }
    public StsFieldBean[] getPropertyFields() { return null; }
	public StsFieldBean[] getDefaultFields() { return null; }
    public boolean anyDependencies() { return true; }
    public Object[] getChildren() { return new Object[0]; }
    public StsObjectPanel getObjectPanel() { return null; }
    public void treeObjectSelected() {}
}
