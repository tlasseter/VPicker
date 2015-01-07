
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;


public class StsVisibilityToolbar extends JPanel implements ActionListener
{
	static private final String wellsLabel = "Wells";
	static private final String surfacesLabel = "Surfaces";
	static private final String pseudosLabel = "Faults/Pseudo Wells";
	static private final String sectionsLabel = "Sections";

	static public final String wellsCommand = "Wells";
	static public final String surfacesCommand = "Surfaces";
	static public final String pseudosCommand = "Faults";
	static public final String sectionsCommand = "Sections";

    private GridBagLayout gridBag = new GridBagLayout();

   	Border border = BorderFactory.createRaisedBevelBorder();
    private JLabel titleLabel  = new JLabel("Display: ");
	private JToggleButton surfaceToggle = new JCheckBox(surfacesLabel);
	private JToggleButton wellToggle = new JCheckBox(wellsLabel);
	private JToggleButton pseudoToggle = new JCheckBox(pseudosLabel);
	private JToggleButton sectionToggle = new JCheckBox(sectionsLabel);
    private JLabel dummyLabel = new JLabel();
	private ActionListener actionListener = null;

    private boolean surfaces;
    private boolean wells;
    private boolean pseudos;
    private boolean sections;

    public StsVisibilityToolbar()
    {
		this(true, true, true, true);
    }

    public StsVisibilityToolbar(boolean surfaces, boolean wells, boolean pseudos, boolean sections)
    {
		this.surfaces = surfaces;
        this.wells = wells;
        this.pseudos = pseudos;
        this.sections = sections;
        try { jbInit(); }
        catch(Exception e)
        {
	       	System.out.println("Exception in : StsVisibilityPanel()\n" + e);
        }
    }

    private void addButton(JToggleButton toggle, String label, boolean selected,
    						String command)
    {
        toggle.setSelected(selected);
        toggle.addActionListener(this);
        toggle.setActionCommand(command);
    }

    private void jbInit() throws Exception
    {
//    	addButton(surfaceToggle, surfacesLabel, surfaces, surfacesCommand);
//    	addButton(wellToggle, wellsLabel, wells, wellsCommand);
    	addButton(pseudoToggle, pseudosLabel, pseudos, pseudosCommand);
    	addButton(sectionToggle, sectionsLabel, sections, sectionsCommand);
        setLayout(gridBag);
        setBorder(border);
        Dimension preferredSize = surfaceToggle.getPreferredSize();
		titleLabel.setFont(new Font("Dialog", 0, 11));
        //preferredSize.height += 8;
        preferredSize.height += 15;
        setPreferredSize(preferredSize);
		add(titleLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4,3,4,3), 3, 0));
//		add(surfaceToggle, new GridBagConstraints2(1, 0, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4,3,4,3), 3, 0));
//		add(wellToggle, new GridBagConstraints2(2, 0, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4,3,4,3), 3, 0));
		add(pseudoToggle, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4,3,4,3), 3, 0));
		add(sectionToggle, new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4,3,4,3), 3, 0));
		add(dummyLabel, new GridBagConstraints(5, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 3, 4, 3), 3, 0));

		setVisible(true);
	}

    public void addActionListener(ActionListener listener)
    {
    	this.actionListener = listener;
    }

    public void setSurfaceToggle(boolean b) { surfaceToggle.setSelected(b); }
    public void setWellsToggle(boolean b) { wellToggle.setSelected(b); }
    public void setPseudosToggle(boolean b) { pseudoToggle.setSelected(b); }
    public void setSectionsToggle(boolean b) { sectionToggle.setSelected(b); }
    public void actionPerformed(ActionEvent e)
    {
//       	System.out.println("VisibilityPanel Command: " + e.getActionCommand());
        JToggleButton toggle = null;
        if( e.getSource() instanceof JToggleButton )
        	toggle = (JToggleButton) e.getSource();

		if( e.getActionCommand().equals(surfacesCommand) )
        {
            surfaces = toggle.isSelected();
	    }
        else if( e.getActionCommand().equals(wellsCommand) )
        {
            wells = toggle.isSelected();
       	}
        else if( e.getActionCommand().equals(pseudosCommand) )
        {
            pseudos = toggle.isSelected();
        }
        else if( e.getActionCommand().equals(sectionsCommand) )
        {
            sections = toggle.isSelected();
        }
        if( actionListener != null )
        	actionListener.actionPerformed(e);
    }

}
