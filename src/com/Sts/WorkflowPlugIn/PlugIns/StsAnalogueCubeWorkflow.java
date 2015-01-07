package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class StsAnalogueCubeWorkflow extends StsWorkflowPlugIn
{
    public StsAnalogueCubeWorkflow()
    {
        workflowName = "Analogue Cube";
        name = "StsAnalogueCubeWorkflow";
        checkName();
        description = new String("Loads, displays, and analyzes 3D volumes using Aspect Energy Intellectual Property.\n" +
                                 " The analysis is referred to as an Analogue Cube and is used to prospect a volume for \n" +
                                 "potentially economic oil & gas reserves\n");

        createComboBoxDescriptors();
    }

    public void createComboBoxDescriptors()
    {
        StsWorkflowPlugIn.ComboBoxDescriptor comboBoxDescriptor;

        String[] volumeClasses = new String[] {"com.Sts.DBTypes.StsSeismicVolume", "com.Sts.DBTypes.StsMathVirtualVolume", "com.Sts.DBTypes.StsBlendedVirtualVolume",
                                           "com.Sts.DBTypes.StsFilterVirtualVolume", "com.Sts.DBTypes.StsRGBAVirtualVolume",
                                           "com.Sts.DBTypes.StsCrossplotVirtualVolume", "com.Sts.DBTypes.StsSensorVirtualVolume" };
        comboBoxDescriptor = constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsSeismicVolume", volumeClasses, "seismic", "noseismic");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
        comboBoxDescriptor = constructComboBoxDescriptor("com.Sts.DBTypes.StsCrossplot", "crossplot", "noCrossplot");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
    }

    public void createWorkflowNodes(StsTreeModel treeModel, StsWorkflowTreeNode root)
    {
        StsNodeBundle nodeBundle;

        treeModel.addMenuNode(newProject);
        treeModel.addMenuNode(openProject);

        StsWorkflowTreeNode workflowRoot = root.addChild("Analog Cube Workflow", null);
        StsWorkflowTreeNode workflowProcessSeismic = workflowRoot.addChild("com.Sts.Actions.Wizards.PostStack.StsPostStack3dWizard", "Process Post-Stack Seismic", "loadSeismic20x20.gif");
        StsWorkflowTreeNode loadData = workflowRoot.addChild("Load Data", null);
        StsWorkflowTreeNode workflowLoadSeismic = loadData.addChild("com.Sts.Actions.Wizards.PostStack3dLoad.StsVolumeWizard", "Post-Stack Seismic", "loadSeismic20x20.gif");
        StsWorkflowTreeNode workflowLoadWells = loadData.addChild("com.Sts.Actions.Wizards.Well.StsWellWizard", "Wells & Logs", "well20x20.gif");
        StsWorkflowTreeNode workflowLoadGrids = loadData.addChild("com.Sts.Actions.Wizards.Surfaces.StsSurfaceWizard", "Surfaces", "importSurfaces20x20.gif");
		StsWorkflowTreeNode loadSpectrums = loadData.addChild("com.Sts.Actions.Wizards.Color.StsPaletteWizard", "Color Palettes", "importPalette20x20.gif");

        StsWorkflowTreeNode workflowSubVolume = workflowRoot.addChild("com.Sts.Actions.Wizards.SubVolume.StsSubVolumeWizard", "Create/Edit SubVolumes", "defineCropVolume20x20.gif");
        StsWorkflowTreeNode workflowHorPick = workflowRoot.addChild("com.Sts.Actions.Wizards.Horpick.StsHorpickWizard", "Pick Surfaces", "horizonPicker20x20.gif");
        StsWorkflowTreeNode workflowDefineHorizons = workflowRoot.addChild("com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Define Horizons", "defineHorizons20x20.gif");
        StsWorkflowTreeNode workflowAnalogueCube = workflowRoot.addChild("com.Sts.Actions.Wizards.AnalogueCube.StsAnalogueCubeWizard", "Create Analog Cube", "newAnalogueCube.gif");

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowProcessSeismic));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, loadSpectrums));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoadWells));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoadGrids));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowLoadSeismic));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowSubVolume));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowLoadSeismic }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowHorPick));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowAnalogueCube));

        logUsageChange();
    }

    public void addOptionalNodes(StsTreeModel treeModel, String[] options)
    {

    }
}
