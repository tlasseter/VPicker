
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to display a list

package com.Sts.UI;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsDoubleListDialog extends JDialog
{
    static public final String PROTOTYPE_STRING = "PROTOTYPE STRING";
    static private final Dimension MIN_SIZE = new Dimension(200, 200);

    private float rowHeight;
    private boolean offset;
    private String button1Text;
    private String button2Text;
    private boolean useCancelBtn;
    private boolean useButton1;
    private boolean useButton2;
    private String label1Text = null;
    protected String[] items = null;
    private String label2Text = null;
    private String[] items2 = null;

	private JPanel mainPanel = new JPanel();
	private BorderLayout paneLayout = new BorderLayout();
	private JPanel buttonPanel = new JPanel();
	private FlowLayout flowLayout1 = new FlowLayout();
    private JLabel label = new JLabel();
    private BorderLayout borderLayout1 = new BorderLayout();
    private JScrollPane listsScrollPane = new JScrollPane();
    private JPanel listsPanel = new JPanel();
    protected JList list2 = new JList();
    private Object itemPrototype;
    private ListCellRenderer cellRenderer;
    private JLabel label2 = new JLabel();
    protected JList list = new JList();
    private JLabel label1 = new JLabel();
    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    protected JButton cancelButton = null;
    protected JButton okayButton = new JButton();
    protected JButton button2 = null;
    protected JButton button1 = null;

	public StsDoubleListDialog(JFrame parent, String title, String text,
            String[] names, float rowHeight)
	{
        this(parent, title, text, null, names, null, null, rowHeight);
    }
	public StsDoubleListDialog(JFrame parent, String title, String text,
            String label1Text, String[] names1, String label2Text, String[] names2,
            float rowHeight)
	{
        this(parent, title, text, label1Text, names1, label2Text, names2, true,
                rowHeight, true, false, false, null, null, PROTOTYPE_STRING, null);
    }
	public StsDoubleListDialog(JFrame parent, String title, String text,
            String label1Text, String[] names1, String label2Text, String[] names2,
            boolean offset, float rowHeight, boolean finishDialog,
            boolean modal, boolean useCancelBtn, String button1Text,
            String button2Text, Object itemPrototype, ListCellRenderer cellRenderer)
	{
        super(parent, title, modal);
        label.setText(text);
        this.label1Text = label1Text;
        this.items = names1;
        this.items2 = names2;
        this.label2Text = label2Text;
        this.rowHeight = rowHeight;
        this.offset = offset;
        this.useCancelBtn = useCancelBtn;
        this.button1Text = button1Text;
        this.button2Text = button2Text;
        this.itemPrototype = itemPrototype;
        this.cellRenderer = cellRenderer;
        useButton1 = (button1Text!=null);
        useButton2 = (button2Text!=null);
		jbInit();
        if (finishDialog)
        {
            setItems(names1, names2);
    		pack();
            adjustSize();
        }
	}

	private void jbInit()
	{
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		mainPanel.setLayout(paneLayout);
		buttonPanel.setLayout(flowLayout1);
        listsPanel.setLayout(gridBagLayout1);
        listsPanel.setBackground(Color.white);
		if (useButton1)
        {
            button1 = new JButton();
            button1.setText(button1Text);
            buttonPanel.add(button1);
        }
		if (useButton2)
        {
            button2 = new JButton();
            button2.setText(button2Text);
            buttonPanel.add(button2);
        }
        okayButton.setText("Okay");
        buttonPanel.add(okayButton);
		if (useCancelBtn)
        {
            cancelButton = new JButton();
            cancelButton.setText("Cancel");
            buttonPanel.add(cancelButton);
        }
        mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        mainPanel.add(label, BorderLayout.NORTH);
        mainPanel.add(listsScrollPane, BorderLayout.CENTER);
        listsScrollPane.getViewport().add(listsPanel, null);
        if (label1Text != null)
        {
            listsPanel.add(label1, new GridBagConstraints(0, 0, 1, 1, 0.5, 0.0,
                    GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
        }
        int height = 0;
        if (itemPrototype != null)
        {
            if (cellRenderer != null) list.setCellRenderer(cellRenderer);
            list.setPrototypeCellValue(itemPrototype);
            height = (int)(list.getFixedCellHeight() * rowHeight);
            list.setFixedCellHeight(height);
        }
        listsPanel.add(list, new GridBagConstraints(0, 1, 1, 1, 0.5, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));
        if (label2Text != null)
        {
            listsPanel.add(label2, new GridBagConstraints(1, 0, 1, 1, 0.5, 0.0,
                    GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
        }
        if (label2Text != null || items2 != null)
        {
            if (cellRenderer != null) list2.setCellRenderer(cellRenderer);
            Insets list2Insets = new Insets(0, 0, 0, 0);
            if (height > 0)
            {
                list2.setFixedCellHeight(height);
                list2Insets.top = height/2;
            }
            listsPanel.add(list2, new GridBagConstraints(1, 1, 1, 1, 0.5, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    list2Insets, 0, 0));
        }
        Border border = BorderFactory.createEtchedBorder();
        if (button1 != null) button1.addActionListener(
            new java.awt.event.ActionListener()
    		{
	    		public void actionPerformed(ActionEvent e)
		    	{
    				button1_actionPerformed(e);
	    		}
		    });
        if (button2 != null) button2.addActionListener(
            new java.awt.event.ActionListener()
    		{
	    		public void actionPerformed(ActionEvent e)
		    	{
    			    	button2_actionPerformed(e);
	    		}
		    });
        okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButton_actionPerformed(e);
			}
		});
        if (cancelButton != null) cancelButton.addActionListener(
            new java.awt.event.ActionListener()
    		{
	    		public void actionPerformed(ActionEvent e)
		    	{
			    	cancelButton_actionPerformed(e);
    			}
	    	});
        if (label1Text != null) label1.setText(label1Text);
        if (label2Text != null) label2.setText(label2Text);
        listsPanel.setBorder(border);
        getContentPane().add(mainPanel);
	}

	protected void okayButton_actionPerformed(ActionEvent e)
	{
		setVisible(false);
    }

	protected void cancelButton_actionPerformed(ActionEvent e)
	{
        items = items2 = null;
        list.clearSelection();
        list2.clearSelection();
		setVisible(false);
	}

	protected void button1_actionPerformed(ActionEvent e)
	{

	}

	protected void button2_actionPerformed(ActionEvent e)
	{

	}

    public void setItems(String[] names)
    {
        setItems(names, null);
    }
    public void setItems(String[] names1, String[] names2)
    {
        items = names1;
        items2 = names2;
        list.setListData(names1);
        list2.setListData(names2);
    }

    public String[] getItems() { return items; }
    public String[] getItems2() { return items2; }
    public JList getList() { return list; }
    public JList getList2() { return list2; }

    protected void adjustSize()
    {
        Dimension dialogSize = getSize();
        if (dialogSize.height < MIN_SIZE.height) dialogSize.height = MIN_SIZE.height;
        if (dialogSize.width < MIN_SIZE.width) dialogSize.width = MIN_SIZE.width;
		Dimension maxSize = Toolkit.getDefaultToolkit().getScreenSize();
        maxSize.height *= 0.75f;
        if (dialogSize.height > maxSize.height) dialogSize.height = maxSize.height;
        this.setSize(dialogSize);
    }

    public static void main(String[] args)
    {
        String[] names1 = { "Red", "Green", "Blue", "Yellow", "Purple", "Black" };
        String[] names2 = { "square", "circle", "triangle", "diamond", "hexagon" };
        StsDoubleListDialog d = new StsDoubleListDialog(new JFrame(), "StsListDialog test",
                "Double List Example", "Colors: ", names1, "Shapes: ", names2, 2);
        d.setModal(true);
        d.setVisible(true);
        System.exit(0);
    }

}

