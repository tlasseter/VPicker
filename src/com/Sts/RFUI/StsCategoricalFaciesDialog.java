
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.Actions.Export.*;
import com.Sts.Help.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsCategoricalFaciesDialog extends JDialog
    implements ActionListener, ListSelectionListener
{
	static public final String OKAY = "Okay";
	static public final String CANCEL = "Cancel";
	static public final String HELP = "Help";

	private StsDoubleListPanel panel = new StsDoubleListPanel();
    private StsCategoricalFacies categoricalFacies;
    private StsExportCategoricalFacies action;
    private Frame frame;
    private String[] faciesNames;

    Border border = BorderFactory.createLoweredBevelBorder();
	JPanel mainPanel = new JPanel();
	GridBagLayout mainLayout = new GridBagLayout();
	JComponent dialogPanel = new JPanel();
	JPanel buttonPanel = new JPanel();
	JButton okayButton = new JButton();
	JButton cancelButton = new JButton();
	JButton helpButton = new JButton();
	FlowLayout buttonLayout = new FlowLayout();

    static final private String HELP_URL = "resframe.html";
    static final private String HELP_TARGET = "top";

	public StsCategoricalFaciesDialog(Frame frame, boolean modal,
            StsCategoricalFacies categoricalFacies,
            StsExportCategoricalFacies action)
	{
        super(frame, "Categorical Facies Assignment Dialog", modal);
        this.categoricalFacies = categoricalFacies;
        this.action = action;
        this.frame = frame;
        faciesNames = categoricalFacies.getFaciesNames();
		jbInit();
 	    pack();
 	}

    private void initPanel()
    {
        panel.setTitle("Assign all lithologies to facies categories:");
        panel.setLeftTitle("Facies Categories:");
        panel.setRightTitle("Lithologies:");
        panel.getLeftButton().setToolTipText("Assign lithologies to a facies category");
        panel.getRightButton().setToolTipText("Remove all lithologies from a facies category");
        panel.setBorder(border);
    	mainPanel.add(panel, new GridBagConstraints(1, 0, 3, 4, 0.75, 1.0
           		    ,GridBagConstraints.EAST, GridBagConstraints.BOTH,
                	new Insets(10, 10, 10, 10), 0, 0));
        panel.addActionListener(this);
        panel.addListSelectionListener(this);

        refresh();
        JList leftList = panel.getLeftList();
        Dimension leftSize = leftList.getPreferredSize();
        JList rightList = panel.getRightList();
        Dimension rightSize = rightList.getPreferredSize();
        int height = Math.max(leftSize.height, rightSize.height);
        leftSize.height = height;
        rightSize.height = height;
        leftList.setSize(leftSize);
        rightList.setSize(rightSize);

        //leftList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }

	private void jbInit()
	{
		mainPanel.setLayout(mainLayout);
        dialogPanel.setBorder(border);
        buttonPanel.setLayout(buttonLayout);
		okayButton.setText(OKAY);
		cancelButton.setText(CANCEL);
		helpButton.setText(HELP);
		okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButton_actionPerformed(e);
			}
		});
		cancelButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				cancelButton_actionPerformed(e);
			}
		});
		helpButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				helpButton_actionPerformed(e);
			}
		});
        HelpManager.setButtonHelp(helpButton, HelpManager.GENERAL, HELP_TARGET, frame);

        initPanel();
		mainPanel.add(buttonPanel, new GridBagConstraints(0, 4, 4, 1, 1.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(10, 0, 10, 0), 0, 0));
		buttonPanel.add(okayButton);
		buttonPanel.add(cancelButton);
		buttonPanel.add(helpButton);
        getContentPane().add(mainPanel);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);

	}

	protected void okayButton_actionPerformed(ActionEvent e)
	{
        if (action != null)
        {
            int[] unassigned = categoricalFacies.getUnassignedLithologies();
            action.setOkay(unassigned == null);
        }
		dispose();
    }

	protected void cancelButton_actionPerformed(ActionEvent e)
	{
        if (action != null) action.setOkay(false);
		dispose();
	}

	protected void helpButton_actionPerformed(ActionEvent e)
	{
//        HelpManager.loadRelativeURL(HELP_URL);
	}

    public void setFaciesList()
    {
        String[] items = new String[faciesNames.length];

        for (int i=0; i<faciesNames.length; i++)
        {
            items[i] = new String(faciesNames[i]);
            int[] lithologies = categoricalFacies.getLithologies(i);
            if (lithologies == null) items[i] += " ( )";
            else
            {
                items[i] += " (";
                for (int j=0; j<lithologies.length; j++)
                {
                    items[i] += (lithologies[j] + " ");
                }
                items[i] += ")";
            }
            int[] selected = panel.getSelectedLeftIndices();
            panel.setLeftItems(items);
            if (selected != null) panel.setSelectedLeftIndices(selected);
        }
    }

    private void setLithologyList()
    {
        String items[] = null;
        int[] unassignedLithologies = categoricalFacies.getUnassignedLithologies();
        if (unassignedLithologies != null)
        {
            items = new String[unassignedLithologies.length];
            for (int i=0; i<unassignedLithologies.length; i++)
            {
                items[i] = new String((new Integer(unassignedLithologies[i])).toString());
            }
        }
        else items = new String[0];
        panel.setRightItems(items);
    }

    public void refresh()
    {
    	setFaciesList();
    	setLithologyList();
        refreshButtons();
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	refreshButtons();
    }

    public void refreshButtons()
    {
		Object[] faciesItems = panel.getSelectedLeftItems();
		Object[] lithologyItems = panel.getSelectedRightItems();
    	int nFaciesItems = (faciesItems == null) ? 0 : faciesItems.length;
    	int nLithologyItems = (lithologyItems == null) ? 0 : lithologyItems.length;
        boolean enableRight = false;
        if (nFaciesItems > 0)
        {
            int facies = StsCategoricalFacies.getFaciesCategory((String)faciesItems[0]);
            if (categoricalFacies.getLithologies(facies) != null) enableRight = true;
        }

        panel.enableLeftButton(nFaciesItems == 1 && nLithologyItems > 0);
        panel.enableRightButton(enableRight);
    }

    public void actionPerformed(ActionEvent e)
    {
    	if (e.getActionCommand().equals("<"))
        {
    		Object[] faciesItems = panel.getSelectedLeftItems();
	    	Object[] lithologyItems = panel.getSelectedRightItems();
      	    int nFaciesItems = (faciesItems == null) ? 0 : faciesItems.length;
        	int nLithologyItems = (lithologyItems == null) ? 0 : lithologyItems.length;
            if (nFaciesItems == 0 || nLithologyItems == 0) return;

            int facies = StsCategoricalFacies.getFaciesCategory((String)faciesItems[0]);
            for (int i=0; i<nLithologyItems; i++)
            {
                String item = (String)lithologyItems[i];
                int lithology = Integer.valueOf(item).intValue();
                categoricalFacies.set(facies, lithology);
             }
             refresh();
        }
        else if (e.getActionCommand().equals(">"))
        {
    		Object[] faciesItems = panel.getSelectedLeftItems();
      	    if (faciesItems == null) return;
            int facies = StsCategoricalFacies.getFaciesCategory((String)faciesItems[0]);
            categoricalFacies.removeFacies(facies);
            refresh();
        }
    }

    public static void main(String[] args)
    {
        StsCategoricalFacies cf = null;
        try { cf = new StsCategoricalFacies(5); }
        catch (StsException e)
        {
            System.out.println("Error creating categorical facies");
            System.exit(0);
        };
        final float[] lithologies = { 1, 2, 3, 4 };
        cf.addLithologies(lithologies);

        /*
        cf.set(1, 1);
        cf.set(2, 2);
        cf.set(3, 3);
        cf.set(4, 4);
        */

        StsCategoricalFaciesDialog d = new StsCategoricalFaciesDialog(null, true, cf, null);
        d.setVisible(true);
        System.exit(0);
    }

}

