
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

public class StsListDialog extends JDialog
{
    private float rowHeight;
    private boolean offset;
    private boolean useCancelBtn;
    protected String[] items = null;

	private JPanel mainPanel = new JPanel();
	private BorderLayout paneLayout = new BorderLayout();
	private JPanel buttonPanel = new JPanel();
	private FlowLayout flowLayout1 = new FlowLayout();
    private JLabel label = new JLabel();
    private BorderLayout borderLayout1 = new BorderLayout();
    private JScrollPane listsScrollPane = new JScrollPane();
    private JPanel listsPanel = new JPanel();
    private Object itemPrototype;
    private ListCellRenderer cellRenderer;
    protected JList list = new JList();
    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    protected JButton cancelButton = null;
    protected JButton okayButton = new JButton();

    static private final Dimension MIN_SIZE = new Dimension(200, 200);
    static public final float DEFAULT_ROW_HEIGHT = 1.1f;
    static public final String PROTOTYPE_STRING = "PROTOTYPE STRING";

	public StsListDialog(JFrame parent, String title, String text, String[] names)
	{
        this(parent, title, text, names, DEFAULT_ROW_HEIGHT, true, false, false,
                StsDoubleListDialog.PROTOTYPE_STRING, null);
    }
	public StsListDialog(JFrame parent, String title, String text,
            String[] names, float rowHeight, boolean finishDialog, boolean modal,
            boolean useCancelBtn, Object itemPrototype, ListCellRenderer cellRenderer)
	{
        super(parent, title, modal);
        label.setText(text);
        this.items = names;
        this.rowHeight = rowHeight;
        this.useCancelBtn = useCancelBtn;
        this.itemPrototype = itemPrototype;
        this.cellRenderer = cellRenderer;
		jbInit();
        if (finishDialog)
        {
            setItems(names);
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

        Border border = BorderFactory.createEtchedBorder();

        okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButtonAction(e);
			}
		});
        if(cancelButton != null)
            cancelButton.addActionListener(new java.awt.event.ActionListener()
    		{
	    		public void actionPerformed(ActionEvent e)
		    	{
			    	cancelButtonAction(e);
    			}
	    	});
        listsPanel.setBorder(border);
        getContentPane().add(mainPanel);
	}


	protected void okayButtonAction(ActionEvent e)
	{
		setVisible(false);
    }

	protected void cancelButtonAction(ActionEvent e)
	{
        items = null;
        list.clearSelection();
		setVisible(false);
	}

    public void setItems(String[] names) { items = names; }

    public String[] getItems() { return items; }
    public JList getList() { return list; }

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
        String[] names = { "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow" };
        StsListDialog d = new StsListDialog(null, "StsListDialog test",
                "list of colors:", names);
        d.setModal(true);
        d.setVisible(true);
        System.exit(0);
    }
}

