package com.Sts.WorkflowPlugIn.PlugIns;

import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class TriconPreStackVelocityWorkflow extends StsWorkflowPlugIn
{
    public TriconPreStackVelocityWorkflow()
    {
        workflowName = "Velocity Analysis";
        name = "TriconPreStackVelocityWorkflow";
        checkName();
        description = new String("Velocity analysis workflow that has been tailored to work more closely with Paradigm Focus." +
        		"Takes Focus DSOUT data as input, then lets user interactively pick stacking velocities" +
        		"using semblance, gather flattening, CVS and VVS panels as guides. Velocities can be output as DISCO Define" +
        		"jobs using HANDVEL cards, which are automatically executed by clicking Finish.");

        createComboBoxDescriptors();
    }

    public void createComboBoxDescriptors()
    {
        StsWorkflowPlugIn.ComboBoxDescriptor comboBoxDescriptor;
        String[] volumeClasses = new String[]
                                            {"com.Sts.DBTypes.StsSeismicVolume",
                                            "com.Sts.DBTypes.StsFilterVirtualVolume" };
        comboBoxDescriptor =
            constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsSeismicVolume",
                    volumeClasses, "poststack", "noseismic");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
        String[] preStackClasses = new String[]
                                              {"com.Sts.DBTypes.StsPreStackLineSet3d",
                                              "com.Sts.DBTypes.StsPreStackLineSet2d"};
        comboBoxDescriptor =
            StsWorkflowPlugIn.constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsPreStackLineSet",
                    preStackClasses, "prestack", "noseismic");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
    }

    public void createWorkflowNodes(StsTreeModel treeModel, StsWorkflowTreeNode workflowRoot)
    {
        treeModel.addMenuNode(newProject);
        treeModel.addMenuNode(openProject);
        
        workflowRoot.setImageIcon("tricon_logo.gif");
        
        StsWorkflowTreeNode importSeismic = workflowRoot.addChild("Import/Reformat Data", "orange_circle_one.gif");
        StsWorkflowTreeNode importFocus3d = importSeismic.addChild("com.tricongeophysics.Actions.Wizards.PreStackFocus3d.StsPreStackFocusWizard3d", "Focus 3D", "focus_logo.gif");
        StsWorkflowTreeNode importFocus2d = importSeismic.addChild("com.tricongeophysics.Actions.Wizards.PreStackFocus2d.StsPreStackFocusWizard2d", "Focus 2D", "focus_logo.gif");
        StsWorkflowTreeNode importSegy3d = importSeismic.addChild("com.Sts.Actions.Wizards.PreStack3d.StsPreStackWizard", "SEG-Y 3D", "loadSeismic20x20.gif");
        StsWorkflowTreeNode importSegy2d = importSeismic.addChild("com.Sts.Actions.Wizards.PreStack2d.StsPreStack2dWizard", "SEG-Y 2D", "2DSeismic20x20.gif");
        
        StsWorkflowTreeNode loadData = workflowRoot.addChild("Load Data", "orange_circle_two.gif");
        StsWorkflowTreeNode workflowLoad3dSeismic = loadData.addChild("com.Sts.Actions.Wizards.PreStack3dLoad.StsPreStackLoadWizard", "3D", "loadSeismic20x20.gif");
        StsWorkflowTreeNode workflowLoad2dSeismic = loadData.addChild("com.Sts.Actions.Wizards.PreStack2dLoad.StsPreStackLoad2dWizard", "2D", "2DSeismic20x20.gif");
        StsWorkflowTreeNode workflowLoadHandVel = loadData.addChild("com.Sts.Actions.Wizards.HandVelocity.StsHandVelocityWizard", "HandVels", "importHandVels20x20.gif");

        StsWorkflowTreeNode buildVelocityModel = workflowRoot.addChild("Build Velocity Model", "orange_circle_three.gif");
        StsWorkflowTreeNode workflowNewVelocity = buildVelocityModel.addChild(
                "com.tricongeophysics.Actions.Wizards.VelocityAnalysis.TriconVelocityAnalysisWizard", 
                "* New *", 
                "Pick Velocities for New Velocity Model",
                "Must Load 2D or 3D dataset before you can pick velocities", 
                "loadSeismic20x20.gif");
        StsWorkflowTreeNode workflowCopyVelocity = buildVelocityModel.addChild("com.Sts.Actions.Wizards.CopyVelocity.StsCopyVelocityWizard", "Copy",   "loadSeismic20x20.gif");
        StsWorkflowTreeNode workflowFilterVelocity = buildVelocityModel.addChild("com.Sts.Actions.Wizards.VirtualVolume.StsFilterVolumeWizard", "Filter (not working)", "defineVirtual20x20.gif");
        
        StsWorkflowTreeNode export = workflowRoot.addChild("Export", "orange_circle_four.gif");
        StsWorkflowTreeNode workflowExport3dHandvel = export.addChild(
                "com.tricongeophysics.Actions.Wizards.PreStackExport.TriconPreStackExportHandvelWizard3d", 
                "3D Handvel DISCO", 
                "Build and Execute Paradigm DISCO Handvel job using current velocity model",
                "Must Build a 3D Velocity Model before you can Export.",
                "importHandVels20x20.gif");
        StsWorkflowTreeNode workflowExport2dHandvel = export.addChild(
                "com.tricongeophysics.Actions.Wizards.PreStackExport.TriconPreStackExportHandvelWizard2d", 
                "2D Handvel DISCO", 
                "Build and Execute Paradigm DISCO Handvel job using current velocity model",
                "Must Build a 2D Velocity Model before you can Export.",
                "importHandVels20x20.gif");
        StsWorkflowTreeNode workflowExportVolume = export.addChild("com.Sts.Actions.Wizards.PreStackExport.StsPreStackExportWizard", "Volumes (not working)", "loadSeismic20x20.gif");
        
        StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, importFocus3d));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, importFocus2d));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, importSegy3d));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, importSegy2d));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoad3dSeismic));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoad2dSeismic));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoadHandVel));
        
        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowLoad3dSeismic, workflowLoad2dSeismic }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowNewVelocity));
        
        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowNewVelocity  }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowCopyVelocity));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowFilterVelocity));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowExportVolume));
        
        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowLoad3dSeismic, workflowNewVelocity  }, StsNodeBundle.ALL_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowExport3dHandvel));
        
        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowLoad2dSeismic, workflowNewVelocity  }, StsNodeBundle.ALL_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowExport2dHandvel));

        logUsageChange();
    }

    public void addOptionalNodes(StsTreeModel treeModel, String[] options)
    {

    }
}
