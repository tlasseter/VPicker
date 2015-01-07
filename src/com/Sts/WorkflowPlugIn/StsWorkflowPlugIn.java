package com.Sts.WorkflowPlugIn;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.*;

import java.util.*;

/** This class describes a plug-in.  A plug-in has a name,
 *  a description, and an optional set of comboBoxDescriptors."
 */
public class StsWorkflowPlugIn
{
    public String name = "com.Sts.MVC.StsMain";
    public String workflowName = "Default";
    public String description = "None";
    public ArrayList comboBoxToolbarDescriptors = new ArrayList();

    StsTreeModel treeModel;
    public StsWorkflowTreeNode root;
    public StsTreeNode rootNode, dataNode, modelNode;
    public StsTreeNode subVolumeNode, virtualVolumeNode, seismicNode, crossplotNode;

    public StsWorkflowTreeNode newProject;
    public StsWorkflowTreeNode openProject;

//    StsWorkflowTreeNode optionalRoot;

    String PROCESS_SEISMIC_GROUP = "Process Seismic Data";
    String LOAD_DATA_GROUP = "Load";
    String DEFINE_EDIT_GROUP = "Define / Edit";
    String ANALYZE_GROUP = "Analyze";
    String OUTPUT_GROUP = "Output";

    public final byte PROCESS_SEISMIC = 0;
    public final byte LOAD_DATA = 1;
    public final byte DEFINE_EDIT = 2;
    public final byte ANALYZE = 3;
    public final byte OUTPUT = 4;

    String[] groupNodeNames = {PROCESS_SEISMIC_GROUP, LOAD_DATA_GROUP, DEFINE_EDIT_GROUP, ANALYZE_GROUP, OUTPUT_GROUP };
    StsWorkflowTreeNode[] groupNodes = new StsWorkflowTreeNode[5];

    public final byte P_PRESTACK2D = 0;
    public final byte P_PRESTACK3D = 1;
    public final byte P_POSTSTACK2D = 2;
    public final byte P_POSTSTACK3D = 3;
    public final byte P_VSP = 4;
    byte[] processDataNodeOrder = { P_PRESTACK2D, P_PRESTACK3D, P_POSTSTACK2D, P_POSTSTACK3D, P_VSP };
    StsWorkflowTreeNode[] processData = { 
    		StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PreStack2d.StsPreStack2dWizard", "PreStack 2D",
    				"Pre-process pre-stack 2D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"2DSeismic20x20.gif"),
        /*
            StsWorkflowTreeNode.constructor("com.tricongeophysics.Actions.Wizards.PreStackFocus2d.StsPreStackFocusWizard2d", "PreStack Focus 2D",
    				"Pre-process Focus pre-stack 2D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"2DSeismic20x20.gif"),
        */
    		StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PreStack3d.StsPreStackWizard", "PreStack 3D",
    				"Pre-process pre-stack 3D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"loadPreSeismic20x20.gif"),
        /*
            StsWorkflowTreeNode.constructor("com.tricongeophysics.Actions.Wizards.PreStackFocus3d.StsPreStackFocusWizard3d", "PreStack Focus 3D",
    				"Pre-process Focus pre-stack 3D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"loadPreSeismic20x20.gif"),
        */
    		StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PostStack2d.StsPostStack2dWizard", "PostStack 2D",
    				"Pre-process post-stack 2D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"2DSeismic20x20.gif"),
    		StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PostStack3d.StsPostStack3dWizard", "PostStack 3D",
    				"Pre-process post-stack 3D datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",
    		"loadSeismic20x20.gif"),
    		StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.Vsp.StsSegyVspWizard", "Vertical Profiles",
    				"Pre-process 2D VSP datasets and optimize them for loading into S2S.",
    				"A project must be created or opened prior to running this workflow step.",    				
    		"loadSeismic20x20.gif")
    };

    public final byte L_PRESTACK2D = 0;
    public final byte L_PRESTACK3D = 1;
    public final byte L_POSTSTACK2D = 2;
    public final byte L_POSTSTACK3D = 3;
    public final byte L_HANDVEL = 4;
    public final byte L_COLORPALETTE = 5;
    public final byte L_WELLSANDLOGS = 6;
    public final byte L_SENSOR = 7;
    public final byte L_VSP = 8;
    public final byte L_CULTURE = 9;
    public final byte L_ANCILLARY = 10;
    public final byte L_SURFACE = 11;
    public final byte L_FAULT_STICKS = 12;
    byte[] loadDataNodeOrder = { L_PRESTACK2D, L_PRESTACK3D, L_POSTSTACK2D, L_POSTSTACK3D, L_HANDVEL, L_COLORPALETTE,
    		L_WELLSANDLOGS, L_SENSOR, L_VSP, L_CULTURE, L_ANCILLARY, L_SURFACE, L_FAULT_STICKS};
    StsWorkflowTreeNode[] loadData = {
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PreStack2dLoad.StsPreStackLoad2dWizard", "PreStack 2D Seismic",
        		"Load previously processed prestack 2d datasets into the project.",
        		"Seismic data must be pre-processed prior to loading. Run the Process Seismic Data->PreStack 2D step first.",
        		"2DSeismic20x20.gif"),
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PreStack3dLoad.StsPreStackLoadWizard", "PreStack 3D Seismic",
        		"Load previously processed prestack 3d datasets into the project.",
        		"Seismic data must be pre-processed prior to loading. Run the Process Seismic Data->PreStack 3D step first.",
        		 "loadSeismic20x20.gif"),
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PostStack2dLoad.StsLine2dWizard", "PostStack 2D Seismic",
        		"Load previously processed poststack 2d datasets into the project.",
        		"Seismic data must be pre-processed prior to loading. Run the Process Seismic Data->PostStack 2D step first.",
        		"2DSeismic20x20.gif"),
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.PostStack3dLoad.StsVolumeWizard", "PostStack 3D Seismic",
        		"Load previously processed poststack 3d datasets into the project.",
        		"Seismic data must be pre-processed prior to loading. Run the Process Seismic Data->PostStack 3D step first.",
        		"loadSeismic20x20.gif"),
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.HandVelocity.StsHandVelocityWizard", "HandVels",
            	"Load handvel formatted velocity profiles into the project. These profiles can then be used to initialize a 2D or 3D pre-stack velocity model by running Define->Velocity Model on the PreStack Velocity Modeling Workflow",
            	"A project must be created or opened prior to loading color palettes.",
        		"importHandVels20x20.gif"),
        	StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.Color.StsPaletteWizard", "Color Palettes",
        		"Load color palettes exported from other applications into the project. Once loaded palettes can be assigned to data objects via the object panel.",
        		"A project must be created or opened prior to loading color palettes.",
        		"importPalette20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.Well.StsWellWizard", "Wells & Logs",
	        	"Load well trajectory, logs, markers, perforations and time-depth functions. All data associated with a single well must be in properly named and formated file(s) to load properly.",
	        	"A project must be created or opened prior to loading wells and logs.",
				"well20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.SensorLoad.StsSensorLoadWizard", "Sensor Data",
		        "Load sensor data. Data needs to be in a column ordered and comma delimited ASCII file with a header row. A minimum of two columns are required containing time and an attribute, position columns are optional.",
		        "A well must be loaded prior to loading sensor data.",
				"timeSensor20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.VspLoad.StsVspLoadWizard", "VSP & Assign Wells",
        		"Load previously processed 2D VSP datasets into the project and assign them to previously loaded wells.",
        		"VSP data must be pre-processed and associated wells must have already been loaded prior to loading VSP data. Run the Process Seismic Data->Vertical Profiles and Load->Wells & Logs workflow steps first.",
				"loadSeismic20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.Culture.StsCultureWizard", "Culture",
	        	"Load culture data in the form of lines, symbols and text. Loaded data is in XML ASCII format, examples are provided in the help.",
	        	"A project must be created or opened prior to loading culture data.",
				"loadCulture20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.AncillaryData.StsAncillaryDataWizard", "Ancillary Data",
		        "Load ancillary data such as documents, spreadsheets, images, movies, etc... Loaded data is simply archived with an icon placed in the 3D view to indicate data is available. The icon can be selected and the file launched.",
		        "A project must be created or opened prior to loading ancillary data.",
				"loadAncillaryData20x20.gif"),
			StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.Surfaces.StsSurfaceWizard", "Surfaces",
			    "Load gridded surface data into the project. Two grid formats are supported and are defined in the help.",
			    "A project must be created or opened prior to loading surface data.",
				"importSurfaces20x20.gif"),
            StsWorkflowTreeNode.constructor("com.Sts.Actions.Wizards.FaultSticks.StsFaultSticksWizard", "Fault Sticks",
	    		"Load ASCII file containing fault sticks"	,
	    		"Must have created an S2S Project.",
	    		"faultSticks20x20.gif")
    };
    StsWorkflowTreeNode[][] nodeGroups = {processData, loadData};
    byte[][] nodeOrder = {processDataNodeOrder, loadDataNodeOrder};

    //StsWorkflowTreeNode processData;
    //StsWorkflowTreeNode loadData;
    StsWorkflowTreeNode analyzeData;
    StsWorkflowTreeNode buildStrat;
    StsWorkflowTreeNode buildFrame;
    StsWorkflowTreeNode buildModel;
/*
    public static int numOptionalSteps = 20;

    public static byte SEGYVOLUME = 0;
    public static byte VIRTUALVOLUME = 1;
    public static byte SUBVOLUME = 2;
    public static byte WELLPLAN = 3;
    public static byte VELOCITYMODEL = 4;
    public static byte LOADPALETTES = 5;
    public static byte LOADSEISMIC = 6;
    public static byte LOADWELLS = 7;
    public static byte LOADSURFACES = 8;
    public static byte ANALOGCUBE = 9;
    public static byte WELLMKRS = 10;
    public static byte HORIZONPICK = 11;
    public static byte MOVIE = 12;
    public static byte EDITTD = 13;
    public static byte DEFINEHORIZONS = 14;
    public static byte DEFINEBOUNDARY = 15;
    public static byte BUILDFRAME = 16;
    public static byte CONSTRUCTMODEL = 17;
    public static byte COMBOVOLUMES = 18;
    public static byte CROSSPLOT = 19;
*/
//    public StsWorkflowTreeNode[] optionalSteps = new StsWorkflowTreeNode[numOptionalSteps];
//    public StsWorkflowTreeNode[] optionalSteps = new StsWorkflowTreeNode[numOptionalSteps];
    public StsNodeBundle projectNodeBundle = null;

    public StsWorkflowPlugIn()
    {
        name = getClass().getName();
        newProject = StsWorkflowTreeNode.constructor("com.Sts.Actions.StsNewModel", "New Project", "newProj20x20.gif");
        openProject = StsWorkflowTreeNode.constructor("com.Sts.Actions.StsOpenModel", "Open Project", null);
        projectNodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
    }

    public String getDescription()
    {
        return description;
    }

    public void logUsageChange()
    {
//        System.out.println("Workflow log - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
//        Main.logUsage("com.Sts.WorkflowPlugIn.PlugIns." + name, "Changed Workflow");
        Main.setLogModule("com.Sts.WorkflowPlugIn.PlugIns." + name, "Changed Workflow");
//        System.out.println("Workflow set to - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
    }

    public ArrayList getComboBoxToolbarDescriptors()
    {
        return comboBoxToolbarDescriptors;
    }

    public static ComboBoxDescriptor constructComboBoxDescriptor(String className, String iconName, String deIconName)
    {
        return new ComboBoxDescriptor(className, iconName, deIconName);
    }

    public static ComboBoxDescriptor constructComboBoxToolbarDescriptor(String parentClass, String[] classNames, String iconName, String deIconName)
    {
        return new ComboBoxDescriptor(parentClass, classNames, iconName, deIconName);
    }

    public static class ComboBoxDescriptor
    {
        public String parentClassName;
        public String[] classNames = null;
        public String selectedIconName;
        public String deselectedIconName;

        public ComboBoxDescriptor(String className, String iconName)
        {
            this(className, iconName, iconName);
        }
        public ComboBoxDescriptor(String className, String selectedIconName, String deselectedIconName)
        {
            this.parentClassName = null;
            this.classNames = new String[] {className};
            this.selectedIconName = selectedIconName;
            this.deselectedIconName = deselectedIconName;
        }
        public ComboBoxDescriptor(String parentClass, String[] classNames, String selectedIconName, String deselectedIconName)
        {
            this.parentClassName = parentClass;
            this.classNames = classNames;
            this.selectedIconName = selectedIconName;
            this.deselectedIconName = deselectedIconName;
        }
    }

    public StsWorkflowTreeNode getNodeNamed(byte group, String nodeName)
    {
    	StsWorkflowTreeNode node = null;
        for(int i=0; i<nodeGroups[group].length; i++)
    	{
    	   if(nodeGroups[group][i].getName().equalsIgnoreCase(nodeName))
    	   {
    	      node = nodeGroups[group][i];
    	      break;
    	   }
    	}
    	return node;
    }

    public StsWorkflowTreeNode getNode(byte group, byte nodeIdx)
    {
    	if(nodeGroups.length < group) return null;
    	if(nodeGroups[group].length < nodeIdx) return null;
    	return nodeGroups[group][nodeIdx];
    }

    public void createGroupAndNodes(StsWorkflowTreeNode root, StsTreeModel treeModel, byte group, byte[] nodes)
    {
    	StsNodeBundle nodeBundle = null;

    	groupNodes[group] = root.addChild(groupNodeNames[group], null);
    	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { newProject, openProject }, StsNodeBundle.ONE_REQUIRED);
    	for(int j=0; j<nodes.length; j++)
    	{
            StsWorkflowTreeNode node = nodeGroups[group][nodes[j]];
            if(node == null) continue;
            groupNodes[group].addChild(node);
    		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, node));
    	}
    }

    public StsWorkflowTreeNode getGroupNode(byte group)
    {
    	return groupNodes[group];
    }

    public void createWorkflowNodes(StsTreeModel treeModel, StsWorkflowTreeNode root)
    {

    }
/*
    public void createOptionalNodes(StsTreeModel treeModel, StsWorkflowTreeNode root)
    {
        this.treeModel = treeModel;
        this.root = root;

        // All Available Workflow Steps Classified into groups.
//        StsWorkflowTreeNode optionalRoot = root.addChild("All Available Steps", null);

        StsWorkflowTreeNode processData = root.addChild("Define Data", null);
            optionalSteps[SEGYVOLUME] = processData.addChild(com.Sts.Actions.Wizards.PostStack.StsPostStackWizard", "Process Seismic", "loadSeismic20x20.gif");
            optionalSteps[VIRTUALVOLUME] = processData.addChild(com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard", "Virtual Volumes", "defineVirtual20x20.gif");
            optionalSteps[COMBOVOLUMES] = processData.addChild(com.Sts.Actions.Wizards.CombinationVolume.StsCombinationVolumeWizard", "Combination Volumes", "defineCombination20x20.gif");
            optionalSteps[SUBVOLUME] = processData.addChild(com.Sts.Actions.Wizards.SubVolume.StsSubVolumeWizard", "Sub-Volumes", "defineCropVolume20x20.gif");

        StsWorkflowTreeNode loadData = root.addChild("Load Data", null);
            optionalSteps[LOADPALETTES] = loadData.addChild(com.Sts.Actions.Wizards.Color.StsPaletteWizard", "Load Palettes", "importPalette20x20.gif");
            optionalSteps[LOADSEISMIC] = loadData.addChild(com.Sts.Actions.Wizards.PostStack3d.StsVolumeWizard", "Load Seismic", "loadSeismic20x20.gif");
            optionalSteps[LOADWELLS] = loadData.addChild(com.Sts.Actions.Wizards.Well.StsWellWizard", "Load Wells & Logs", "well20x20.gif");
            optionalSteps[LOADSURFACES] = loadData.addChild(com.Sts.Actions.Wizards.Surfaces.StsSurfaceWizard", "Load Surfaces", "importSurfaces20x20.gif");

        StsWorkflowTreeNode analyzeData = root.addChild("Analyze Data", null);
            optionalSteps[ANALOGCUBE] = analyzeData.addChild(com.Sts.Actions.Wizards.AnalogueCube.StsAnalogueCubeWizard", "Create Analog Cube", "newAnalogueCube.gif");
            optionalSteps[CROSSPLOT] = analyzeData.addChild(com.Sts.Actions.Wizards.Crossplot.StsCrossplotWizard", "Select/Create a Cross Plot", "newCrossplot.gif");
            optionalSteps[HORIZONPICK] = analyzeData.addChild(com.Sts.Actions.Wizards.Horpick.StsHorpickWizard", "Pick Surfaces", "horizonPicker20x20.gif");
            optionalSteps[MOVIE] = analyzeData.addChild(com.Sts.Actions.Wizards.Movie.StsMovieWizard", "Animations", "createMovie20x20.gif");

        StsWorkflowTreeNode buildStrat = root.addChild("Build Strat", null);
            optionalSteps[DEFINEHORIZONS] = buildStrat.addChild(com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Define Horizons", "defineHorizons20x20.gif");
            optionalSteps[DEFINEBOUNDARY] = buildStrat.addChild(com.Sts.Actions.Build.StsBuildBoundary", "Define Boundary", "boundary20x20.gif");

        optionalSteps[BUILDFRAME] = root.addChild(com.Sts.Actions.Build.StsBuildFrame", "Build Frame", "buildFrame20x20.gif");
        optionalSteps[BUILDFRAME].isOptional(true);

        StsWorkflowTreeNode buildModel = root.addChild("Build Model", null);
            optionalSteps[CONSTRUCTMODEL] = buildModel.addChild(com.Sts.Actions.Wizards.Model.StsModelWizard", "Complete model", "buildModel20x20.gif");

        addConnection(SEGYVOLUME, projectNodeBundle);
        addConnection(LOADPALETTES, projectNodeBundle);
        addConnection(LOADSEISMIC, projectNodeBundle);
        addConnection(LOADWELLS, projectNodeBundle);
        addConnection(LOADSURFACES, projectNodeBundle);
        addConnection(MOVIE, projectNodeBundle);

        addConnection(VIRTUALVOLUME, optionalSteps[LOADSEISMIC]);
        addConnection(COMBOVOLUMES, optionalSteps[LOADSEISMIC]);
        addConnection(SUBVOLUME, optionalSteps[LOADSEISMIC]);
        addConnection(ANALOGCUBE, optionalSteps[LOADSEISMIC]);
        addConnection(HORIZONPICK, optionalSteps[LOADSEISMIC]);
        addConnection(CROSSPLOT, optionalSteps[LOADSEISMIC]);

        StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { optionalSteps[HORIZONPICK], optionalSteps[LOADSURFACES] }, StsNodeBundle.ONE_REQUIRED);
        addConnection(DEFINEHORIZONS, nodeBundle);
        addConnection(DEFINEBOUNDARY, optionalSteps[DEFINEHORIZONS]);
        addConnection(BUILDFRAME, optionalSteps[DEFINEBOUNDARY]);
        addConnection(CONSTRUCTMODEL, optionalSteps[BUILDFRAME]);

        loadData.setSelectedNode();

        logUsageChange();
    }

    public void resetConnections(byte step, StsNodeBundle nodeBundle)
    {
        treeModel.removeNodeConnections(optionalSteps[step]);
        addConnection(step, nodeBundle);
    }

    public void addConnection(byte step, StsWorkflowTreeNode node)
    {
        treeModel.addNodeConnection(new StsNodeConnection(node, optionalSteps[step]));
    }

    public void addConnection(byte step, StsNodeBundle nodeBundle)
    {
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, optionalSteps[step]));
    }

    public void addOptionalNodes(StsTreeModel treeModel, String[] options)
    {
    }
*/
    public boolean checkName()
    {
        String className = getClass().getName();
        if(className.indexOf(name) > 0) return true;

        StsException.systemError(className + " constructor incorrect. name: " + name +
            " must be a substring of className" + className);
        return false;
    }

    public StsTreeNode getSubVolumeNode()
    {
        return subVolumeNode;
    }

    public StsTreeNode[] getObjectsThatSupportSubVolumes()
    {
        return new StsTreeNode[] { seismicNode, virtualVolumeNode, crossplotNode };
    }

    public StsTreeNode getNodeForType(Object objectType)
    {
        // Traverse the root tree to find the node where the user object equals the supplied object
        StsTreeNode node = findTreeNode(rootNode, objectType);
        return node;
    }

    public StsTreeNode findTreeNode(StsTreeNode _root, Object objectType)
    {
        StsTreeNode node = null;
        Enumeration children = dataNode.children();
        if (children != null)
        {
            while (children.hasMoreElements())
            {
                node = findTreeNode(_root, (StsTreeNode) children.nextElement());
                if(node != null) return node;
            }
        }
        // Traverse the tree to find the node where the user object equals the supplied object
        return null;
    }

    /** This is a standard objectsPanel.  Override this method to create customized versions. */
	public void createObjectsPanel(StsObjectTree objectTree, StsModel model)
	{
        final StsObjectTree _objectTree = objectTree;
		final StsModel _model = model;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                runCreateObjectsPanel(_objectTree, _model);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    public void addAdditionalToolbars(StsWin3dBase win3d)
    {
        ; // Override in specific workflows as required.
    }

    public void addTimeActionToolbar(StsWin3d win3d)
    {
        if(win3d == null) return;
        win3d.checkAddTimeActionToolbar();
    }

    public StsTreeNode checkAddStaticNode(StsModel model, StsTreeNode node, StsTreeObjectI userObject, String label, boolean required)
    {
        return node.addStaticNode(userObject, label);
        /*
        if(required)
            return node.addDynamicNode(userObject, label);
        // If not required but data of type is already loaded, create node
        else
        {
            if(model.hasObjectsOfType(userObject.getClass(), StsParameters.NONE))
                return
            else
                return null;
        }
        */
    }
    public StsTreeNode checkAddDynamicNode(StsModel model, StsTreeNode node, StsSubType userObject, String label, boolean required)
    {
        return node.addDynamicNode(userObject, label);
    }

    public StsTreeNode checkAddDynamicNode(StsModel model, StsTreeNode node, StsClass userObject, String label, boolean required)
    {
        if(node == null)
            return null;

        if(required)
            return node.addDynamicNode(userObject, label);
        // If not required but data of type is already loaded, create node
        else
        {
            if(userObject.hasObjects())
                return node.addDynamicNode(userObject, label);
            else
                return null;
        }

        // Add the workflow wizards to the workflow help
    }

    protected boolean runCreateObjectsPanel(StsObjectTree objectTree, StsModel model)
    {
        try
        {
            if(objectTree == null) return false;
            rootNode = objectTree.createRootNode(model.getProject(), "Project");
            rootNode.addStaticNode(model.getProject().getCropVolume(), "Crop Project");

            dataNode = rootNode.addStaticNode("Data");
            modelNode = rootNode.addStaticNode("Model");

            seismicNode = dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume"), "3D Volumes");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers");
			dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPatchVolume"), "Curvature Patches");
//            dataNode.addDynamicNode(model.getCreateStsClassDebug("com.Sts.DBTypes.StsPreStackLineSet2d"), "PreStack3d 2D Volumes");
			StsVspClass vspClass = (StsVspClass)model.getCreateStsClass("com.Sts.DBTypes.StsVsp");
            StsTreeNode vspNode = dataNode.addDynamicNode(vspClass, "VSP");
			StsSubType[] subTypes = vspClass.getSubTypes();
			for(int n = 0; n < subTypes.length; n++)
				vspNode.addDynamicNode(subTypes[n], subTypes[n].getName());
/*
            for(int i=0; i< StsVsp.nodeStrings.length; i++)
            {
                String classStg;
				if (i == 0)
                classStg = "com.Sts.DBTypes.StsVsp";
			    else
			    classStg = "com.Sts.DBTypes.StsVsp" + StsVsp.nodeStrings[i];
                vspNode.addDynamicNode(model.getCreateStsClass(classStg), StsVsp.nodeStrings[i]);
            }
 */
            //virtualVolumeNode = dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes");
            virtualVolumeNode = dataNode.addStaticNode(model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMathVirtualVolume"), "Math");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsBlendedVirtualVolume"), "Blended");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsCrossplotVirtualVolume"), "Crossplot");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsRGBAVirtualVolume"), "RGBA");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsFilterVirtualVolume"), "Filter");
            virtualVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSensorVirtualVolume"), "Sensor");

            crossplotNode = dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsCrossplot"), "Crossplot");

            subVolumeNode = dataNode.addStaticNode(model.getCreateStsClass("com.Sts.DBTypes.StsSubVolume"), "SubVolumes");
            subVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsDualSurfaceSubVolume"), "Dual Surface");
            subVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsBoxSetSubVolume"), "Box Set");
            subVolumeNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWellSubVolume"), "Well Set");

            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types");

            StsTreeNode sensorNode = dataNode.addStaticNode(model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensor");
            sensorNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static");
            sensorNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic");

			dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsTank"), "Storage Tanks");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPump"), "Pumps");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitored");

            StsTreeNode alarmNode = dataNode.addStaticNode(model.getCreateStsClass("com.Sts.DBTypes.StsAlarm"), "Alarm");
            alarmNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSurfaceAlarm"), "Surface");
            alarmNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWellAlarm"), "Well");
            alarmNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsValueAlarm"), "Value");

            StsAncillaryDataClass adClass = (StsAncillaryDataClass)model.getCreateStsClass("com.Sts.DBTypes.StsAncillaryData");
            StsTreeNode adNode = dataNode.addDynamicNode(adClass, "Ancillary Data");

            subTypes = adClass.getSubTypes();
            for(int n = 0; n < subTypes.length; n++)
				adNode.addDynamicNode(subTypes[n], subTypes[n].getName());

            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWeighPoint"), "WayPoints");

            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsModelSurface"), "Horizons");
			modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d"), "PreStack3d 3D Velocity Model");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel2d"), "PreStack3d 2D Velocity Model");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVelocityModel"), "Velocity Model");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsEditTdSet"), "Well TD Edits");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsLine"), "Boundary Lines");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsFaultLine"), "Fault Lines");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsZone"), "Units");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSection"), "Sections");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsFractureSet"), "Fracture Sets");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsBlock"), "Blocks");
            modelNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsBuiltModel"), "Sim Models");

            objectTree.finalizeTreeModel();
            return true;

        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowPlug.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
/*
	public void createTreeModel(StsModel model)
    {
        rootNode = StsTreeNode.staticNode(model.getProject(), "Project");

        StsTreeNode subVolumeNode = StsTreeNode.staticNode(model.getProject().getCropVolume(), "Crop PostStack3d");
        rootNode.add(subVolumeNode);
        StsTreeNode dataNode = StsTreeNode.staticNode(null, "Data");
        rootNode.add(dataNode);
        StsTreeNode modelNode = StsTreeNode.staticNode(null, "Model");
        rootNode.add(modelNode);

        // Actual Objects from files or databases
        StsTreeNode seismicVolumesNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSeismicVolume.class), "Seismic Volumes");
        dataNode.add(seismicVolumesNode);
        StsTreeNode surfacesNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSurface.class), "Surfaces");
        dataNode.add(surfacesNode);
        StsTreeNode wellsNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsWell.class), "Wells");
        dataNode.add(wellsNode);
        StsTreeNode spectrumsNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSpectrum.class), "Spectrums");
        dataNode.add(spectrumsNode);

        // Constructed Objects from other objects
        StsTreeNode crossplotNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsCrossplot.class), "Cross Plots");
        dataNode.add(crossplotNode);
        StsTreeNode virtualVolumesNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsVirtualVolume.class), "Virtual Volumes");
        dataNode.add(virtualVolumesNode);

        StsTreeNode dataSubVolumeNode = StsTreeNode.staticNode(null, "SubVolumes");
        dataNode.add(dataSubVolumeNode);
        StsTreeNode polyhedronNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSubVolume.class), "Polyhedron");
        dataSubVolumeNode.add(polyhedronNode);
        StsTreeNode singleSurfaceNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSubVolume.class), "Single Surface");
        dataSubVolumeNode.add(singleSurfaceNode);
        StsTreeNode dualSurfaceNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsSubVolume.class), "Dual Surface");
        dataSubVolumeNode.add(dualSurfaceNode);

        // Constructed Objects from the environment
        StsTreeNode moviesNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsMovie.class), "Movies");
        dataNode.add(moviesNode);
        StsTreeNode weighPointsNode = StsTreeNode.dynamicNode(model.getCreateStsClass(StsWeighPoint.class), "Weigh Points");
        dataNode.add(weighPointsNode);



        treeModel = new DefaultTreeModel(rootNode);
		rootNode.explore();
        tree.setModel(treeModel);
        tree.setSelectionRow(0);
	}
*/
   public void testWorkflow(StsModel model) {}

}