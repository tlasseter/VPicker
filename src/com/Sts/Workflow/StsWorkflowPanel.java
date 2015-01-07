package com.Sts.Workflow;

import com.Sts.Actions.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;
import java.util.*;

/**
 * Title:        WorkFlow
 * Description:  Develop a simple workflow from JTree
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author Tom Lasseter
 * @version 1.0
 */

/** contains the tree which defines a sequence of workflow steps.
 *  Each step has an associated action Class which is executed when the step is selected.
 *  Each step has a state and a colored light which indicates whether the step
 *  CANNOT_START (red), CAN_START (yellow), STARTED (green) meaning it can executed again,
 *  or ENDED (blue) meaning it cannot be executed again.
 *  @see StsWorkflowTreeNode
 *  @see StsNodeConnection
 *  @see StsNodeBundle
 */
public class StsWorkflowPanel extends StsJPanel
{
    /** the controller for all workflow actions */
	StsActionManager actionManager;
    /** the graphical tree displaying the workflow steps */
	JTree tree;
    /** the underlying model for the tree */
    StsTreeModel treeModel;
    /** the workflow combobox list */
    StsComboBoxFieldBean workflowsBean = null;
    String[] plugInNames = null;
    String selectedPlugin = null;
    boolean selectedPluginStatus = true;
    
    /** Contains descriptions of the workflows and workflow steps. */
    JScrollPane detailPane = null;
    JScrollPane prerequisitePane = null;
    JTabbedPane tabbedPane = null;
    
    final byte DETAIL = 0;
    final byte PREREQUISITE = 1;
    
    /**Construct the frame*/
//    public StsWorkflowPanel(StsActionManager actionManager, StsWorkflowPlugIn selectedPlugIn)
    public StsWorkflowPanel(StsModel model, StsActionManager actionManager, StsWorkflowPlugIn selectedPlugIn, String[] options)
	{
        try
		{
            this.actionManager = actionManager;
            // Create tree
//            createTree(model, selectedPlugIn);
            // Build Current Workflow Tree
            createTree(model, selectedPlugIn, options);
/*
            treeModel = StsTreeModel.construct(this, model, selectedPlugIn, options);
		    tree = new JTree(treeModel);
            StsWorkflowTreeNode.setTreeModel(treeModel);
            initializeActionStatus();

			// Put tree in a scrollable pane
			JScrollPane treePane = new JScrollPane(tree);
            treePane.setPreferredSize(new Dimension(150, 400));
            JScrollPane detailPane = new JScrollPane();

	        JSplitPane split = new JSplitPane();
            split.setOrientation(0);
            split.setOneTouchExpandable(true);

			// use a layout that will stretch tree to panel size
			setLayout(new BorderLayout());

            add(split, BorderLayout.CENTER);
		    split.add(treePane, JSplitPane.TOP);
		    split.add(detailPane, JSplitPane.BOTTOM);

			// make customized cell renderer and editor
		    StsWorkflowTreeNodeRenderer renderer = new StsWorkflowTreeNodeRenderer();
			renderer.setClosedIcon(null);
			renderer.setOpenIcon(null);
			renderer.setLeafIcon(null);
		    StsWorkflowTreeNodeEditor editor = new StsWorkflowTreeNodeEditor(actionManager);

		    tree.setEditable(true);
		    tree.setCellRenderer(renderer);
		    tree.setCellEditor(new ImmediateEditor(tree, renderer, editor));

			// adjust tree parameters
			tree.setRowHeight(25);
          */
        }
        catch(Exception e)
		{
            e.printStackTrace();
        }
    }
    
    private void setWorkflowSelectComboBoxItems(StsModel model)
    {
        try
        {
            String[] workflowPlugInNames = model.workflowPlugInNames;
            if(workflowPlugInNames == null) return;
            int nPlugIns = workflowPlugInNames.length;

            plugInNames = new String[nPlugIns];
            for(int n = 0; n < nPlugIns; n++)
            	plugInNames[n] = workflowPlugInNames[n].substring(workflowPlugInNames[n].lastIndexOf(".")+1, 
            			workflowPlugInNames[n].length());
        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowPanel.setWorkflowSelectComboBoxItems() failed.",
                e, StsException.WARNING);
        }
    }
    
    public void clearTree()
    {
        if(tree != null)  tree = null;
        if(treeModel != null) treeModel = null;
    }

//    public void createTree(StsModel model, StsWorkflowPlugIn selectedPlugIn)
    public void createTree(StsModel model, StsWorkflowPlugIn selectedPlugIn, String[] options)
    {
//        treeModel = StsTreeModel.construct(this, model, selectedPlugIn);
        // Add Workflow Combo Box
        if(selectedPlugIn == null) return;
        
        setWorkflowSelectComboBoxItems(model);
        workflowsBean = new StsComboBoxFieldBean(this, "workflowPlugin", "", plugInNames);
        
        treeModel = StsTreeModel.construct(this, model, selectedPlugIn, options);
        tree = new JTree(treeModel);
        StsWorkflowTreeNode.setTreeModel(treeModel);
        initializeActionStatus();

        // Put tree in a scrollable pane
        JScrollPane treePane = new JScrollPane(tree);
        treePane.setPreferredSize(new Dimension(150, 300));
        detailPane = new JScrollPane();
        prerequisitePane = new JScrollPane();
        
        JSplitPane split = new JSplitPane();
        split.setOrientation(0);
        split.setOneTouchExpandable(true);
        
        tabbedPane = new JTabbedPane();
        tabbedPane.addTab("Details", detailPane);
        tabbedPane.addTab("PreRequisites", prerequisitePane);

        // use a layout that will stretch tree to panel size
        //setLayout(new BorderLayout());

        //add(split, BorderLayout.CENTER);
        gbc.weighty = 0.0;
        gbc.fill = gbc.HORIZONTAL;
        addEndRow(workflowsBean);
        gbc.weighty = 1.0;
        gbc.gridwidth = 2;
        gbc.fill = gbc.BOTH;
        addEndRow(split);
                
        updateDetails(actionManager.getModel().getWorkflowPlugIn().getDescription());
        
        split.add(treePane, JSplitPane.TOP);
        split.add(tabbedPane, JSplitPane.BOTTOM);
        split.setResizeWeight(0.5);
        // make customized cell renderer and editor
        StsWorkflowTreeNodeRenderer renderer = new StsWorkflowTreeNodeRenderer();
        renderer.setClosedIcon(null);
        renderer.setOpenIcon(null);
        renderer.setLeafIcon(null);
        StsWorkflowTreeNodeEditor editor = new StsWorkflowTreeNodeEditor(this);

        tree.setEditable(true);
        tree.setCellRenderer(renderer);
        tree.setCellEditor(new ImmediateEditor(tree, renderer, editor));

        // Adjust tree parameters
        tree.setRowHeight(25);
        
        int index = model.getWorkflowPlugInIndex();
        if (index >= 0)
        {
        	workflowsBean.setSelectedIndex(index);
        }
        if(model.win3d != null)
            selectedPlugIn.createObjectsPanel(model.win3d.objectTreePanel, model);
        model.refreshObjectPanel();
    }
    
    public void updateDetails(String text, Color background, Font font)
    {    	
    	JTextArea typeDescription = new JTextArea(text);
        typeDescription.setLineWrap(true);
        typeDescription.setWrapStyleWord(true);
        typeDescription.setBackground(background);
        typeDescription.setFont(font);        
        detailPane.getViewport().add(typeDescription);
        validate();
    }
    
    public void updatePrerequisites(String text, Color background, Font font)
    {    	
    	JTextArea typeDescription = new JTextArea(text);
        typeDescription.setLineWrap(true);
        typeDescription.setWrapStyleWord(true);
        typeDescription.setBackground(background);
        typeDescription.setFont(font);        
        prerequisitePane.getViewport().add(typeDescription);
        validate();
    }
    
    public void updateDetails(String text)
    {
        if(selectedPluginStatus)
    	    updateDetails(text, Color.WHITE, new java.awt.Font("SansSerif", 0, 11));
        else
            updateDetails(text + " --- NOT A LICENSED WORKFLOW.", Color.RED, new java.awt.Font("SansSerif", 0, 11));
    }

    public void updatePrerequisites(String text, boolean failed)
    {    	
    	if(failed)
    	{
    		updatePrerequisites(text, Color.RED, new java.awt.Font("SansSerif", 1, 11));
    		tabbedPane.setSelectedIndex(PREREQUISITE);
    	}
    	else
    		updatePrerequisites(text, Color.WHITE, new java.awt.Font("SansSerif", 0, 11));
    }
    
    public void setWorkflowPlugin(String plugin)
    {
        try
        {
        	selectedPlugin = plugin;
        	StsModel model = actionManager.getModel();
            ClassLoader classLoader = getClass().getClassLoader();
            if((model.workflowPlugInNames == null) || (workflowsBean == null))
            	return;
            int workflowIndex = workflowsBean.getSelectedIndex();
            selectedPluginStatus = model.workflowPlugInStatus[workflowIndex].booleanValue();
            Class plugInClass = classLoader.loadClass(model.workflowPlugInNames[workflowIndex]);
            Constructor constructor = plugInClass.getConstructor(new Class[0]);
            StsWorkflowPlugIn selectedPlugIn = (StsWorkflowPlugIn)constructor.newInstance(new Class[0]);
            model.setWorkflowPlugIn(selectedPlugIn);
            updateDetails(actionManager.getModel().getWorkflowPlugIn().getDescription());
        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowPanel.plugInSelected() failed.",
                                         e, StsException.WARNING);
        }
    } 
    public boolean getWorkflowPluginStatus()
    {
    	return selectedPluginStatus;
    }
    public String getWorkflowPlugin()
    {
    	return selectedPlugin;
    }

/*
    public String[] getOptionalSteps()
    {
        return options;
    }

    public void addOptionalStep(String step)
    {
        options = (String[])StsMath.arrayAddElement(options, step);
    }

    public void removeOptionalStep(String step)
    {
        options = (String[])StsMath.arrayDeleteElement(options,step);
    }
*/

    /** initializeActionStatus
     *
     *  Actions are leafs on the workflow tree display.  The color of the light indicates the actionStatus.
     *
     *  CANNOT_START (red) - required input not available
     *  CAN_START (yellow) - required input available, but step has not been executed
     *  STARTED (green) - step has been executed one or more times
     *  ENDED (blue) - step has been terminated because: 1) it can be executed only once,
     *                 2) a next step using this one has been executed thus locking it from further execution
     *
     *  pass 1: For steps in the CANNOT_START state, set them to CAN_START if input is now available or
     *          there are no input requirements.
     *          A single pass is sufficient since a new CAN_START does not propagate downstream.
     *  pass 2: For Tree display, do a bottom-up breadth-first enumeration. For each parent (not a leaf),
     *          set the actionStatus based on that of the children.  Find the maximum and minimum states
     *          for the children. If all in the same state, set parent to that state; otherwise set
     *          parent to the minimum state.
     *
     */
    public void initializeActionStatus()
    {
        StsWorkflowTreeNode root, node;
        Enumeration enumeration;

        // set status for menu nodes
        if(treeModel == null) return;
        StsWorkflowTreeNode[] menuNodes = treeModel.menuNodes;
        for(int n = 0; n < menuNodes.length; n++)
            menuNodes[n].initializeActionStatus();

        // set status for each node in the tree
        root = (StsWorkflowTreeNode)treeModel.getRoot();
        enumeration = root.breadthFirstEnumeration();
        while(enumeration.hasMoreElements())
        {
            node = (StsWorkflowTreeNode)enumeration.nextElement();
            node.initializeActionStatus();
        }

        // now adjust state of leaf nodes based on states of other nodes
        // if changed, update parent status (non-leaf) which will recursively
        // adjust up the chain
        enumeration = root.breadthFirstEnumeration();
        while(enumeration.hasMoreElements())
        {
            node = (StsWorkflowTreeNode)enumeration.nextElement();
            if(node.updateActionStatus()) node.updateParentStatus();
        }

        // set workflow selection button to on
        root.setActionStatus(root.CAN_START);
    }
/*
    private void addMouseListener()
    {
        MouseListener ml = new MouseAdapter()
        {
            public void mousePressed(MouseEvent e)
            {
                int selRow = tree.getRowForLocation(e.getX(), e.getY());
                TreePath path = tree.getPathForRow(selRow);
                if(path == null) return;
                StsWorkflowTreeNode treeNode = (StsWorkflowTreeNode)path.getLastPathComponent();

                Class actionClass = treeNode.getActionClass();
                if(actionClass == null || actionManager == null) return;
                int actionStatus = treeNode.getActionStatus();
                if(actionStatus == StsWorkflowTreeNode.CANNOT_START)
                {
                    StsMessageFiles.logMessage("Action: " + treeNode.name + " can't be executed: input not available.");
                }
                else if(actionStatus == StsWorkflowTreeNode.ENDED)
                {
                    StsMessageFiles.logMessage("Action: " + treeNode.name + " can't be executed: already executed or bypassed.");
                }
                else
                    actionManager.startAction(actionClass);
/*
//                 TreePath selPath = tree.getPathForLocation(e.getX(), e.getY());
                 if(selRow != -1)
                 {
                     if(e.getClickCount() == 1)
                     {
//                         mySingleClick(selRow, selPath);
                     }
                     else if(e.getClickCount() == 2)
                     {
//                         myDoubleClick(selRow, selPath);
                     }
                 }
*/
/*
            }
        };
        tree.addMouseListener(ml);
    }
*/
    public void adjustWorkflowPanelState(String actionClassName, int actionStatus)
    {
        treeModel.adjustWorkflowPanelState(actionClassName, actionStatus);
    }

    /**Main method*/
    public static void main(String[] args)
	{
        int width = 300;
        int height = 200;
		String title = "Workflow Test";

        try
		{
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            StsWorkflowPanel workflowPanel = new StsWorkflowPanel(null, null, null,null);
			JFrame frame = new JFrame(title);
			frame.setSize(width, height);
			Container contentPane = frame.getContentPane();
		    contentPane.add(workflowPanel, BorderLayout.CENTER);

			Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
			Dimension frameSize = frame.getSize();
			if (frameSize.height > screenSize.height)
			{
				frameSize.height = screenSize.height;
			}
			if (frameSize.width > screenSize.width)
			{
				frameSize.width = screenSize.width;
			}
			frame.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
            frame.addWindowListener
			(
				new WindowAdapter()
				{
			        public void windowClosing(WindowEvent e) {System.exit(0); }
				}
			);
			frame.setVisible(true);
        }
        catch(Exception e)
		{
            e.printStackTrace();
        }
    }
}

/** inner class renders the workflow step button with light, icon, and label */
class StsWorkflowTreeNodeRenderer extends DefaultTreeCellRenderer
{
    JPanel panel;
    JLabel iconLabel;
    JButton actionButton;
    JCheckBox endCheckBox;
    JLabel statusLabel;
    JPanel leftButtonsPanel;

    static final int CANNOT_START = StsModel.CANNOT_START;
    static final int CAN_START = StsModel.CAN_START;
    static final int STARTED = StsModel.STARTED;
    static final int ENDED = StsModel.ENDED;

	static ImageIcon[] statusIcons = new ImageIcon[4];
	static
	{
		statusIcons[0] = StsIcon.createIcon("Workflow/RedSmall20x20.gif");    // CANNOT_START
		statusIcons[1] = StsIcon.createIcon("Workflow/YellowSmall20x20.gif"); // CAN_START
		statusIcons[2] = StsIcon.createIcon("Workflow/GreenSmall20x20.gif");  // STARTED
		statusIcons[3] = StsIcon.createIcon("Workflow/BlueSmall20x20.gif");  // ENDED
	}

    ImageIcon statusIcon = statusIcons[0];

    public StsWorkflowTreeNodeRenderer()
	{
		panel = new JPanel();
		// Uncomment this line and comment the line after to reimplement the button look on the panel - JKF
        //panel.setBorder(BorderFactory.createRaisedBevelBorder());
		  panel.setOpaque(false);
//        panel.setMinimumSize(new Dimension(190, 30));
//        panel.setPreferredSize(new Dimension(190, 30));
        panel.setLayout(new BorderLayout());

		iconLabel = new JLabel();
//        iconLabel.setAlignmentX((float) 0.5);
//        iconLabel.setMaximumSize(new Dimension(24, 24));
//        iconLabel.setMinimumSize(new Dimension(24, 24));
        iconLabel.setOpaque(false);
//        iconLabel.setPreferredSize(new Dimension(24, 24));
//		iconLabel.setIcon(createIcon("well20x20.gif"));
//		iconLabel.setIcon(StsIcon.createIcon("well20x20.gif"));

		actionButton = new JButton();
        actionButton.setFont(new Font("Dialog", Font.PLAIN, 11));
        actionButton.setMinimumSize(new Dimension(140, 26));
        actionButton.setPreferredSize(new Dimension(140, 26));
//        actionButton.setLayout(borderLayout);
        actionButton.setText("Action Button Label");

		leftButtonsPanel = new JPanel();
        leftButtonsPanel.setLayout(new BorderLayout());
        leftButtonsPanel.setOpaque(false);
//        leftButtonsPanel.setMaximumSize(new Dimension(40, 26));
//        rightButtonsPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 2, 2));
//        rightButtonsPanel.setLayout(borderLayout);
//        rightButtonsPanel.setOpaque(false);
//        rightButtonsPanel.setMaximumSize(new Dimension(40, 26));

		endCheckBox = new JCheckBox();
		endCheckBox.setSelected(false);
//        endCheckBox.setOpaque(false);
        endCheckBox.setMaximumSize(new Dimension(24, 24));
//        endCheckBox.setMinimumSize(new Dimension(24, 24));
//        endCheckBox.setPreferredSize(new Dimension(24, 24));

		statusLabel = new JLabel();
//        statusLabel.setMaximumSize(new Dimension(24, 24));
//        statusLabel.setMinimumSize(new Dimension(24, 24));
		//statusLabel.setOpaque(true);
		statusLabel.setOpaque(false);
//        statusLabel.setPreferredSize(new Dimension(24, 24));
		statusLabel.setIcon(statusIcon);
//	    statusLabel.setOpaque(false);

        leftButtonsPanel.add(statusLabel, BorderLayout.WEST);
        leftButtonsPanel.add(iconLabel, BorderLayout.CENTER);
        leftButtonsPanel.add(endCheckBox, BorderLayout.EAST);
        panel.add(leftButtonsPanel, BorderLayout.WEST);
//        panel.add(statusLabel, BorderLayout.WEST);
        panel.add(actionButton, BorderLayout.CENTER);
//        panel.add(iconLabel, BorderLayout.EAST);
//       panel.add(rightButtonsPanel, BorderLayout.EAST);
    }

    public void setTooltip(String tooltip)
    {
        if(actionButton != null)
        {
            actionButton.setToolTipText(tooltip);
        }
    }
	public Component getTreeCellRendererComponent(JTree tree, Object value,
						boolean selected, boolean expanded, boolean leaf,
						int row, boolean hasFocus)
	{
		StsWorkflowTreeNode node = (StsWorkflowTreeNode)value;
//        panel.setBackground(node.getColor());
		iconLabel.setIcon(node.getImageIcon());
        actionButton.setText(node.getName());
        actionButton.setHorizontalAlignment(SwingConstants.LEFT);

        int actionStatus = node.getActionStatus();

        if(node.isOptional)
        {
            endCheckBox.setVisible(true);
            endCheckBox.setSelected(actionStatus == ENDED);
        }
        else
            endCheckBox.setVisible(false);
//		endCheckBox.setSelected(node.isActionEnded());
//		boolean isSelected = node.isNodeSelected();
/*
		if(node.isNodeSelected())
			panel.setBackground(Color.blue);
		else
			panel.setBackground(Color.lightGray);
*/
        statusIcon = statusIcons[node.getActionStatus()];
        statusLabel.setIcon(statusIcon);

//        panel.validate();

		super.getTreeCellRendererComponent(tree, value, selected, expanded, leaf, row, hasFocus);
		return panel;
	}

    public Rectangle getCheckBoxBounds(int x, int y)
    {
        Rectangle box = endCheckBox.getBounds();
        box.translate(x, y);
        return box;
    }
/*
	static private ImageIcon createIcon(String name)
	{
		try
		{
		    Class thisClass = Class.forName("com.Sts.Workflow.StsWorkflowTreeNodeRenderer");
			String pathName = new String("Icons/" + name);
            java.net.URL url = thisClass.getResource(pathName);
            if(url == null)
            {
                System.out.println("StsWorkflowPanel.StsWorkflowTreeNodeRenderer.createIcon() error.  Couldn't find image: " + "Icons/" + name);
                return null;
            }
            else
            {
			    Image i = Toolkit.getDefaultToolkit().createImage(url);
			    return new ImageIcon(i);
            }
		}
		catch (Exception ex)
        {
            StsException.systemError("StsWorkflowPanel.createIcon() failed.\n" +
                "Icon file: " + "Icons/" + name);
            return null;
        }
    }
*/
}

class StsWorkflowTreeNodeEditorRenderer extends StsWorkflowTreeNodeRenderer
{
	public Component getTreeCellRendererComponent(
						JTree tree, Object value,
						boolean selected, boolean expanded,
						boolean leaf, int row,
						boolean hasFocus)
	{
		Component c = super.getTreeCellRendererComponent(tree,
							value, selected, expanded,
							leaf, row, hasFocus);
//		setIcon(null);
		return c;
	}

	public JButton getActionButton() { return actionButton; }
	public JCheckBox getEndCheckBox() { return endCheckBox; }
}

class StsAbstractCellEditor implements TableCellEditor, TreeCellEditor
{
    protected EventListenerList listenerList = new EventListenerList();
	protected Object value;
    protected ChangeEvent changeEvent = null;
    protected int clickCountToStart = 1;

	public Object getCellEditorValue()
	{
		return value;
	}

	public void setCellEditorValue(Object value)
	{
		this.value = value;
	}

    public void setClickCountToStart(int count)
	{
		clickCountToStart = count;
    }

    public int getClickCountToStart()
	{
		return clickCountToStart;
    }

    // must return true if cell is to be selectable (mousePressed)
    public boolean isCellEditable(EventObject anEvent)
	{
	/*
		if (anEvent instanceof MouseEvent)
		{
	    	if (((MouseEvent)anEvent).getClickCount() < clickCountToStart)
			    return false;
		}
	*/
		return true;
    }

    public boolean shouldSelectCell(EventObject anEvent)
	{
		if (this.isCellEditable(anEvent))
		{
	    	if (anEvent == null || ((MouseEvent)anEvent).getClickCount() >= clickCountToStart)
				return true;
		}
		return false;
	}

    public boolean stopCellEditing()
	{
		fireEditingStopped();
		return true;
    }
    public void cancelCellEditing()
	{
		fireEditingCanceled();
    }
    public void addCellEditorListener(CellEditorListener l)
	{
		listenerList.add(CellEditorListener.class, l);
    }

    public void removeCellEditorListener(CellEditorListener l)
	{
		listenerList.remove(CellEditorListener.class, l);
    }

    public Component getTreeCellEditorComponent ( JTree tree, Object value,
				boolean isSelected, boolean expanded, boolean leaf, int row)
	{
		return null;
    }

    public Component getTableCellEditorComponent( JTable table, Object value,
						 boolean isSelected, int row, int column)
	{
		return null;
    }

    protected void fireEditingStopped()
	{
		Object[] listeners = listenerList.getListenerList();

		for (int i = listeners.length-2; i>=0; i-=2)
		{
		    if (listeners[i]==CellEditorListener.class)
			{
				if (changeEvent == null)
		    		changeEvent = new ChangeEvent(this);

				((CellEditorListener)listeners[i+1]).editingStopped(changeEvent);
	    	}
		}
    }

    protected void fireEditingCanceled()
	{
		Object[] listeners = listenerList.getListenerList();

		for (int i = listeners.length-2; i>=0; i-=2) {
		    if (listeners[i]==CellEditorListener.class)
			{
				if (changeEvent == null)
		    		changeEvent = new ChangeEvent(this);

			((CellEditorListener) listeners[i+1]).editingCanceled(changeEvent);
		    }
		}
    }
}

class StsWorkflowTreeNodeEditor extends StsAbstractCellEditor
{
	StsWorkflowTreeNodeEditorRenderer renderer;
	StsWorkflowTreeNode lastEditedNode;
	StsActionManager actionManager;
	StsWorkflowPanel workflowPanel;
	JButton actionButton;
	JCheckBox endCheckBox;

	public StsWorkflowTreeNodeEditor(StsWorkflowPanel panel)
	{
		this.actionManager = panel.actionManager;
		workflowPanel = panel;
		renderer = new StsWorkflowTreeNodeEditorRenderer();
		actionButton = renderer.getActionButton();
		actionButton.addMouseListener(
				new MouseListener()
		{
			public void mousePressed(MouseEvent e)
			{
				if (e.getButton() == MouseEvent.BUTTON1)
				{
					showDescription();
					showPrerequisites();
					startAction();
				}
			}

			public void mouseReleased(MouseEvent e)
			{
			}

			public void mouseEntered(MouseEvent e)
			{
                showDescription();
			}

			public void mouseExited(MouseEvent e)
			{
			}

			public void mouseClicked(MouseEvent e)
			{
			}
		});
		/*actionButton.addActionListener
		(
			new ActionListener()
			{
			    public void actionPerformed(ActionEvent e)
			    {
				    startAction();
			    }
		    }
		);*/

        endCheckBox = renderer.getEndCheckBox();
        if(endCheckBox == null) return;
		endCheckBox.addActionListener
		(
			new ActionListener()
			{
			    public void actionPerformed(ActionEvent e)
				{
                    if(lastEditedNode.isCheckBoxSelectable())
				        lastEditedNode.setActionStatusEnded(endCheckBox.isSelected());
                    else
                    {
                        StsMessageFiles.logMessage("Action: " + lastEditedNode.name + " can't be executed: input not available.");
                        Toolkit.getDefaultToolkit().beep();
                    }
				}
			}
		);
	}
	public void showDescription()
	{
		workflowPanel.updateDetails(lastEditedNode.getNodeDescription());
	}
	public void showPrerequisites()
	{
		workflowPanel.updatePrerequisites(lastEditedNode.getNodeRequirements(), false);		
	}	
	public void showPrerequisites(String reqs)
	{
		workflowPanel.updatePrerequisites(reqs, true);
	}	
	public void startAction()
	{
        //System.out.println("Starting action: " + lastEditedNode.getActionClass().getName() + " Status= " + lastEditedNode.getActionStatus());
        if((actionManager == null) || (lastEditedNode == null)) return;

        // Check to see if the action is already running.
        if(actionManager.getCurrentAction() != null)
        {
            if(lastEditedNode.getActionClass().getName().contains(actionManager.getCurrentAction().getName()))
            {
                StsMessageFiles.logMessage(actionManager.getCurrentAction().getName() + " is already running. If you double-clicked on workflow step, only one mouse click is required.");
                return;
            }
        }
        Class actionClass = lastEditedNode.getActionClass();
        int actionType = lastEditedNode.getActionType();

		if(actionClass == null)
        {
            StsMessageFiles.logMessage("Action: " + lastEditedNode.name + " is a summary workstep. Pick an actual work-step.");
            
            // Expand the group
            actionManager.getModel().win3d.getWorkflowPanel().tree.expandPath(actionManager.getModel().win3d.getWorkflowPanel().tree.getSelectionPath());
            Toolkit.getDefaultToolkit().beep();
            return;
        }
        int actionStatus = lastEditedNode.getActionStatus();
        if(workflowPanel.selectedPluginStatus == false)
        {
            StsMessageFiles.logMessage("Action: NOT LICENSED FOR THIS WORKFLOW or WORKFLOW STEP.");
            showPrerequisites("NOT LICENSED FOR THIS WORKFLOW or WORKFLOW STEP. Contact S2S Systems at sales@s2ssystems.com to get access.");
            Toolkit.getDefaultToolkit().beep();
        }
        else if(actionStatus == StsWorkflowTreeNode.CANNOT_START)
        {
            StsMessageFiles.logMessage("Action: " + lastEditedNode.name + " can't be executed: input not available.");
            showPrerequisites(lastEditedNode.getNodeRequirements());
            Toolkit.getDefaultToolkit().beep();
        }
        else if(actionStatus == StsWorkflowTreeNode.ENDED)
        {
            StsAction action = actionManager.newWorkflowAction(actionClass, lastEditedNode.getName());
            if(actionManager.checkStartAction(action) == null)
                StsMessageFiles.logMessage("Action: " + lastEditedNode.name + " not execute on user decision.");
        }
        else
        {
            if(actionType == StsWorkflowTreeNode.ACTIVE)
            {
                StsAction action = actionManager.newWorkflowAction(actionClass, lastEditedNode.getName());
                if(actionManager.checkStartAction(action) == null)
                {
                	if(action.getReasonForFailure() != null)
                	{
                		showPrerequisites(action.getReasonForFailure());
                		Toolkit.getDefaultToolkit().beep();
                	}
                	return;
                }
            }
            else
            {
                try
                {
                    Object[] arguments = new Object[] { actionManager };
                    Class[] argTypes = new Class[] { StsActionManager.class };
                    Constructor constructor = actionClass.getDeclaredConstructor(argTypes);
                    StsWizard passiveWizard = (StsWizard)constructor.newInstance(arguments);

                    if(passiveWizard instanceof Runnable)
                    {
                        Thread thread = new Thread((Runnable)passiveWizard);
                        thread.start();
                        System.out.println("PassiveWizard thread " + thread.getName() + " started.");
                    }
                    else
                        passiveWizard.start();
                }
                catch(Exception e)
                {
                    StsMessageFiles.logMessage("Passive action: " + lastEditedNode.name + " can't be executed." + e);
                    Toolkit.getDefaultToolkit().beep();
                }
            }
        }
	}

	public Component getTreeCellEditorComponent(JTree tree, Object value,
						boolean selected, boolean expanded, boolean leaf, int row)
	{
		lastEditedNode = (StsWorkflowTreeNode)value;
		return renderer.getTreeCellRendererComponent(tree, value, selected,
							expanded, leaf, row, true); // hasFocus ignored
	}

    public boolean isCellEditable(EventObject anEvent)
	{
        return true;
    /*
        if(lastEditedNode == null) return true;
        int actionStatus = lastEditedNode.getActionStatus();
        return actionStatus >= StsWorkflowTreeNode.CAN_START;
    */
    }

    public boolean shouldSelectCell(EventObject anEvent)
	{
		if (this.isCellEditable(anEvent))
		{
	    	if (anEvent == null || ((MouseEvent)anEvent).getClickCount() >= clickCountToStart)
				return true;
		}
		return false;
	}
/*
    // must return true if cell is to be selectable (mousePressed)
    public boolean isCellEditable(EventObject anEvent)
	{
        return lastEditedNode.actionStatus >= StsWorkflowTreeNode.STARTED;
    }

    public boolean shouldSelectCell(EventObject anEvent)
	{
        return isCellEditable(anEvent);
	}
*/

	public Object getCellEditorValue()
	{
		return lastEditedNode.getUserObject();
	}
}

class ImmediateEditor extends DefaultTreeCellEditor
{
	private StsWorkflowTreeNodeRenderer renderer;

	public ImmediateEditor(JTree tree,
							StsWorkflowTreeNodeRenderer renderer,
							StsWorkflowTreeNodeEditor editor)
	{
		super(tree, renderer, editor);
		this.renderer = renderer;
	}

    protected boolean canEditImmediately(EventObject event)
    {
        if(inCheckBoxHitRegion((MouseEvent)event)) return true;

        if((event instanceof MouseEvent) && SwingUtilities.isLeftMouseButton((MouseEvent)event))
        {
            MouseEvent  me = (MouseEvent)event;
            return ((me.getClickCount() > 0) && inHitRegion(me.getX(), me.getY()));
        }
        return (event == null);
    }

	protected boolean XcanEditImmediately(EventObject e)
	{
		if(e instanceof MouseEvent)
		{
			MouseEvent me = (MouseEvent)e;
			TreePath path = tree.getPathForLocation(me.getX(), me.getY());
			StsWorkflowTreeNode node = (StsWorkflowTreeNode) path.getLastPathComponent();
//            DefaultTreeModel treeModel = (DefaultTreeModel)tree.getModel();
//            treeModel.nodeChanged(node);
//            node.setSelectedNode();
            return true;
        /*
			if(node.isNodeSelected())
				return true;
			else
			{
				node.setSelectedNode();
				return false;
			}
        */
		}
		return false;
	}

    public boolean shouldSelectCell(EventObject event)
    {
        if(inCheckBoxHitRegion((MouseEvent)event)) return false;
//        return true;
//        return false;
	    return realEditor.shouldSelectCell(event);
    }

    public boolean XshouldSelectCell(EventObject e)
	{
		if(e instanceof MouseEvent)
		{
			MouseEvent me = (MouseEvent)e;
			TreePath path = tree.getPathForLocation(me.getX(), me.getY());
			StsWorkflowTreeNode node = (StsWorkflowTreeNode) path.getLastPathComponent();
			return node.isLeaf();
		}
		else
		    return false;
	}

    public boolean inCheckBoxHitRegion(MouseEvent me)
    {
        TreePath path = tree.getPathForLocation(me.getX(), me.getY());
        StsWorkflowTreeNode node = (StsWorkflowTreeNode) path.getLastPathComponent();
        if(!node.isCheckBoxSelectable()) return false;
        Rectangle bounds = tree.getRowBounds(lastRow);
		if(bounds == null || renderer == null) return false;
        Rectangle checkBoxBounds = renderer.getCheckBoxBounds(bounds.x, bounds.y);
        Point point = me.getPoint();
        return checkBoxBounds.contains(point);
    }  

}

// A prefix of nodes with an associated constraint used to determine if input requirements for
// a node are satisfied.
/*
class NodeBundle
{
    StsWorkflowTreeNode[] nodes;
    int constraint;

    static final int NONE_REQUIRED = 0;
    static final int ONE_REQUIRED = 1;
    static final int ALL_REQUIRED = 2;

    NodeBundle(StsWorkflowTreeNode[] nodes, int constraint)
    {
        this.nodes = nodes;
        if(nodes == null) nodes = new StsWorkflowTreeNode[0];
        this.constraint = constraint;
    }

    NodeBundle(StsWorkflowTreeNode node)
    {
        if(node == null) nodes = new StsWorkflowTreeNode[0];
        else             nodes = new StsWorkflowTreeNode[] { node };
        constraint = ALL_REQUIRED;
    }

    protected void addOutputNode(StsWorkflowTreeNode outNode)
    {
        int nNodes = nodes.length;
        for(int n = 0; n < nNodes; n++)
            nodes[n].addOutputNode(outNode);
    }

    protected boolean constraintSatisfied()
    {
        StsWorkflowTreeNode node;

        if(constraint == NONE_REQUIRED) return true;

        boolean satisfied = true;
        int nNodes = nodes.length;
        int nSatisfied = 0;
        for(int n = 0; n < nNodes; n++)
            if(nodes[n].getActionStatus() > StsWorkflowTreeNode.CAN_START) nSatisfied++;

        if(constraint == ONE_REQUIRED) return nSatisfied > 0;
        else return nSatisfied == nNodes;
    }

    protected boolean adjustedState(String actionClassName, int actionStatus)
    {
        int nNodes = nodes.length;
        for(int n = 0; n < nNodes; n++)
            if(nodes[n].actionClassName == actionClassName)
            {
                nodes[n].setActionStatus(actionStatus);
                return true;
            }
        return false;
    }

    protected void setActionStatus(int actionStatus)
    {
        for(int n = 0; n < nodes.length; n++)
            nodes[n].setActionStatus(actionStatus);
    }
}
*/
/*
//  A WorkFlowTreeNode connected to this WorkFlowTreeNode

class NodeConnection
{
    NodeBundle inNodes;
    StsWorkflowTreeNode outNode;

    NodeConnection(StsWorkflowTreeNode inNode, StsWorkflowTreeNode outNode)
    {
        inNodes = new NodeBundle(inNode);
        this.outNode = outNode;
        inNodes.addOutputNode(outNode);
        outNode.addInputNodeBundle(inNodes);
    }

    NodeConnection(NodeBundle inNodes, StsWorkflowTreeNode outNode)
    {
        this.inNodes = inNodes;
        this.outNode = outNode;
        inNodes.addOutputNode(outNode);
        outNode.addInputNodeBundle(inNodes);
    }

    protected boolean adjustWorkflowPanelState(String actionClassName, int actionStatus)
    {
        boolean changed = false;
        if(inNodes.adjustedState(actionClassName, actionStatus))
        {
            if(outNode.updateActionStatus()) changed = true;
        }
        else if(outNode.adjustedState(actionClassName, actionStatus))
        {
            if(outNode.actionStatus > StsWorkflowTreeNode.STARTED)
            {
                inNodes.setActionStatus(StsWorkflowTreeNode.ENDED);
                changed = true;
            }

        }
        return changed;
    }
}
*/