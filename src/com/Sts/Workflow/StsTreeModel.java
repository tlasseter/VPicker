package com.Sts.Workflow;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;
import com.Sts.Help.HelpManager;

import javax.swing.tree.*;

public class StsTreeModel extends DefaultTreeModel
{
    public StsModel model;
    public StsNodeConnection[] nodeConnections = null;
    public StsWorkflowTreeNode[] menuNodes = null;

    public StsTreeModel(StsWorkflowPanel workflowPanel, StsModel model, StsWorkflowTreeNode root, StsWorkflowPlugIn selectedPlugIn)
    {
        super(root);
        this.model = model;
        StsWorkflowTreeNode.setModel(model);
        if(selectedPlugIn == null) return;
        selectedPlugIn.createWorkflowNodes(this, root);
/*
        ArrayList plugIns = model.getPlugIns();
        for(int n = 0; n < plugIns.size(); n++)
        {
            StsWorkflowPlugIn plugIn = (StsWorkflowPlugIn)plugIns.get(n);
            plugIn.createWorkflowNodes(this, root);
        }
*/
    }

    public StsTreeModel(StsWorkflowPanel workflowPanel, StsModel model, StsWorkflowTreeNode root, StsWorkflowPlugIn selectedPlugIn, String[] options)
    {
        super(root);
        this.model = model;
        StsWorkflowTreeNode.setModel(model);
        if(selectedPlugIn == null) return;
        selectedPlugIn.createWorkflowNodes(this, root);

        //        selectedPlugIn.addOptionalNodes(this, options);
/*
        ArrayList plugIns = model.getPlugIns();
        for(int n = 0; n < plugIns.size(); n++)
        {
            StsWorkflowPlugIn plugIn = (StsWorkflowPlugIn)plugIns.get(n);
            plugIn.createWorkflowNodes(this, root);
        }
*/
    }
/*
	static public StsTreeModel construct(StsWorkflowPanel workflowPanel, StsModel model, StsWorkflowPlugIn selectedPlugIn)
	{
        StsWorkflowTreeNode root = new StsWorkflowTreeNode("Select Workflow...", null);
        StsTreeModel treeModel = new StsTreeModel(workflowPanel, model, root, selectedPlugIn);
        return treeModel;
    }
*/
    static public StsTreeModel construct(StsWorkflowPanel workflowPanel, StsModel model, StsWorkflowPlugIn selectedPlugIn, String[] options)
    {
    	// Combobox that has available workflows.
    	
        //StsWorkflowTreeNode root = new StsWorkflowTreeNode(com.Sts.Actions.Wizards.Workflow.StsWorkflowWizard.class, "Select Workflow... ", null);
        if(selectedPlugIn == null) return null;
        StsWorkflowTreeNode root = new StsWorkflowTreeNode(selectedPlugIn.workflowName + " Workflow", null);
        StsTreeModel treeModel = new StsTreeModel(workflowPanel, model, root, selectedPlugIn, options);
        treeModel.initializeHelp(model);
        return treeModel;
    }

    public void initializeHelp(StsModel model)
    {
        try
        {

            if(model.win3d == null) return;

            HelpManager manager = HelpManager.getHelpManager(model.win3d);
            manager.initWorkflowHelpSets(model.win3d.getWorkflowHelpButton());

            for(int i=0; i<nodeConnections.length; i++)
            {
                String hsName = nodeConnections[i].outNode.actionClassName;
                hsName = "Wizards/" + hsName.substring(hsName.lastIndexOf(".")+1) + ".hs";
                if(manager.addWorkflowHelpSet(hsName))
                    StsMessageFiles.infoMessage("Added helpset for wizard: " + hsName);
                else
                    StsMessageFiles.errorMessage("Unable to load helpset for wizard: " + hsName);
            }
        }       
        catch(Exception ex)
        {
            StsMessageFiles.errorMessage("Unable to initialize the workflow help.");
            ex.printStackTrace();
        }
    }

//    public StsWorkflowTreeNode getRoot() { return (StsWorkflowTreeNode)root; }
    public StsNodeConnection[] getNodeConnections() { return nodeConnections; }
    public void setNodeConnections(StsNodeConnection[] nodeConnections) { this.nodeConnections = nodeConnections; }
    public StsWorkflowTreeNode[] getMenuNodes() { return menuNodes; }
    public void setMenuNodes(StsWorkflowTreeNode[] menuNodes) { this.menuNodes = menuNodes; }

    public void addMenuNode(StsWorkflowTreeNode menuNode)
    {
        menuNodes = (StsWorkflowTreeNode[])StsMath.arrayAddElement(menuNodes, menuNode);
    }

    public void addNodeConnection(StsNodeConnection nodeConnection)
    {
        nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, nodeConnection);
    }

    public void removeNodeConnections(StsWorkflowTreeNode node)
    {
        if(nodeConnections == null) return;
        int nStsNodeConnections = nodeConnections.length;
        for(int n = 0; n < nStsNodeConnections; n++)
        {
            if(nodeConnections[n].outNode == node)
            {
                nodeConnections = (StsNodeConnection[]) StsMath.arrayDeleteElement(nodeConnections, n);
                break;
            }
        }

    }

/*
    private void createNodes(StsWorkflowTreeNode root)
    {
        StsNodeBundle nodeBundle;

        StsWorkflowTreeNode newProject = new StsWorkflowTreeNode(StsNewModel.class, "New Project", "newProj20x20.gif");
        StsWorkflowTreeNode openProject = new StsWorkflowTreeNode(StsOpenModel.class, "Open Project", null);
        menuNodes = new StsWorkflowTreeNode[] { newProject, openProject };

//      StsWorkflowTreeNode defineProject = root.addChild("Define project", null);
//	StsWorkflowTreeNode reOpenProject = defineProject.addChild(StsSelectProjectDir.class, "Reopen Project", "reopen20x20.gif");
//      StsWorkflowTreeNode newProject = defineProject.addChild(StsNewModel.class, "New Project", "newProj20x20.gif");
//      StsWorkflowTreeNode loadProject = defineProject.addChild(StsOpenModel.class, "Load project", "openDB20x20.gif");

//        StsWorkflowTreeNode defineData = root.addChild("Define Data", null);
//            StsWorkflowTreeNode defineColor = defineData.addChild(StsColorWizard.class, "Color Palette", "loadSeismic20x20.gif");
//
//            nodeBundle = new NodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, NodeBundle.ONE_REQUIRED);
//            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, defineColor));

        StsWorkflowTreeNode processData = root.addChild("Process Data", null);
            StsWorkflowTreeNode processSeismic = processData.addChild(StsPostStackWizard.class, "Process Seismic", "loadSeismic20x20.gif");

            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, processSeismic));

        StsWorkflowTreeNode loadData = root.addChild("Load Data", null);
            StsWorkflowTreeNode loadSpectrums = loadData.addChild(StsPaletteWizard.class, "Load Palettes", "importPalette20x20.gif");
            StsWorkflowTreeNode loadSeismic = loadData.addChild(StsVolumeWizard.class, "Load Seismic", "loadSeismic20x20.gif");
            StsWorkflowTreeNode loadWells = loadData.addChild(StsWellWizard.class, "Load Wells", "well20x20.gif");
            StsWorkflowTreeNode loadGrids = loadData.addChild(StsSurfaceWizard.class, "Load Surfaces", "importSurfaces20x20.gif");
//		StsWorkflowTreeNode loadFaults = loadData.addChild(StsImportFaultCuts.class, "Load Faults", "importFaults20x20.gif");

            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, loadSpectrums));
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, loadSeismic));
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, loadWells));
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, loadGrids));

        StsWorkflowTreeNode analyzeData = root.addChild("Analyze Data", null);
            StsWorkflowTreeNode crossplot = analyzeData.addChild(StsCrossplotWizard.class, "Select/Create a Cross Plot", "newCrossplot.gif");
//            StsWorkflowTreeNode crossplot = analyzeData.addChild(StsCrossplotSetup.class, "Select/Create a Cross Plot", "activateXplot20x20.gif");
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { loadSeismic }, StsNodeBundle.ONE_REQUIRED);
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, crossplot));


//        StsWorkflowTreeNode outputData = root.addChild("Output", null);
//            StsWorkflowTreeNode createMovies = analyzeData.addChild(StsMovieWizard.class, "Create Movies", "createMovie20x20.gif");
//
//            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
//            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, createMovies));

    if(!Main.isCoreLabs)
    {
        StsWorkflowTreeNode buildStrat = root.addChild("Build Strat", null);
            StsWorkflowTreeNode defineHorizons = buildStrat.addChild(StsHorizonsWizard.class, "Define Horizons", "defineHorizons20x20.gif");
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] {loadGrids}, StsNodeBundle.ONE_REQUIRED);
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(nodeBundle, defineHorizons));

            StsWorkflowTreeNode buildBoundary = buildStrat.addChild(StsBuildBoundary.class, "Define Boundary", "boundary20x20.gif");
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(defineHorizons, buildBoundary));


        StsWorkflowTreeNode buildFrame = root.addChild(StsBuildFrame.class, "Build Frame", "buildFrame20x20.gif");
            buildFrame.isOptional(true);
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(buildBoundary, buildFrame));

        StsWorkflowTreeNode buildModel = root.addChild("Build Model", null);

            StsWorkflowTreeNode completeModel = buildModel.addChild(StsModelWizard.class, "Complete model", "buildModel20x20.gif");
            nodeConnections = (StsNodeConnection[])StsMath.arrayAddElement(nodeConnections, new StsNodeConnection(buildFrame, completeModel));
    }

        loadData.setSelectedNode();
	}
*/
    public void adjustWorkflowPanelState(String actionClassName, int actionStatus)
    {
        StsNodeConnection nodeConnection;
        StsWorkflowTreeNode inNode, outNode;

        if(nodeConnections == null) return;
        int nStsNodeConnections = nodeConnections.length;
        boolean changed = false;
        for(int n = 0; n < nStsNodeConnections; n++)
            if(nodeConnections[n].adjustWorkflowPanelState(actionClassName, actionStatus)) changed = true;

//        if(changed) tree.validate();
//        if(changed) tree.updateUI();
//        if(changed) SwingUtilities.updateComponentTreeUI(tree);
/*
        if(changed)
        {
            tree.invalidate();
            tree.validate();
            tree.repaint();
        }
*/
//        if(changed) treeModel.reload();
//        if(changed) tree.treeDidChange();
//        if(changed) tree.repaint();
    }
}
