package com.Sts.UI;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsListSelectionDialog extends JDialog
{
    static public final String PROTOTYPE_STRING = "PROTOTYPE STRING";
    static private final Dimension MIN_SIZE = new Dimension(200, 200);
    // private int[] selectedIndices = null;
    private float rowHeight;
    private String button1Text;
    private String button2Text;
    private boolean useCancelBtn;
    private boolean useButton1;
    private String label1Text = null;
    protected Object[] items = null;
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
    private JLabel label1 = new JLabel();
    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    protected JButton cancelButton = null;
    protected JButton okayButton = new JButton();
    protected JButton button1 = null;
    static public final float DEFAULT_ROW_HEIGHT = 1.1f;

    public StsListSelectionDialog(Frame parent, String title, String text, Object[] names)
    {
        this(parent, title, text, null, names, 1.0f);
    }

    public StsListSelectionDialog(Frame parent, String title, String text, Object[] names, float rowHeight)
    {
        this(parent, title, text, null, names, rowHeight);
    }

    public StsListSelectionDialog(Frame parent, String title, String text, String label1Text, Object[] names1, float rowHeight)
    {
        this(parent, title, text, label1Text, names1, rowHeight, true, false, false, null, PROTOTYPE_STRING, null);
    }

    private void doFinishDialog()
	{
		setItems(this.items);
		pack();
		adjustSize();
    }

	private void finishDialog()
	{
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { doFinishDialog(); }});
	}

    public StsListSelectionDialog(Frame parent, String title, String text, String label1Text, Object[] names1,
        float rowHeight, boolean finishDialog, boolean modal, boolean useCancelBtn, String button1Text, Object itemPrototype,
        ListCellRenderer cellRenderer)
        {
            super(parent, title, modal);
            label.setText(text);
            this.label1Text = label1Text;
            this.items = names1;
            this.rowHeight = rowHeight;
            this.useCancelBtn = useCancelBtn;
            this.button1Text = button1Text;
            this.itemPrototype = itemPrototype;
            this.cellRenderer = cellRenderer;
            useButton1 = (button1Text != null);
            jbInit();
            if (finishDialog)
            {
				finishDialog();
            }
    }

    static public Object getSingleSelectFromListDialog(JFrame parent, String title, String text, Object[] objects)
    {
        StsListSelectionDialog selectDialog = new StsListSelectionDialog(parent, title, text, objects);
        selectDialog.setModal(true);
        selectDialog.setSingleSelect();
        StsToolkit.centerComponentOnScreen(selectDialog);
        selectDialog.setVisible(true);
        int[] selectedIndices = selectDialog.getSelectedIndices();
        if(selectedIndices == null)
            return null;
        int nSelected = selectedIndices.length;
        if(nSelected == 0) return null;
        return objects[selectedIndices[0]];
    }

    static public Object[] getMultiSelectFromListDialog(Frame parent, String title, String text, Object[] objects)
    {
        return getMultiSelectFromListDialog(parent, parent, title, text, objects);
    }

    static public Object[] getMultiSelectFromListDialog(Frame parent, Component centerComponent, String title, String text, Object[] objects)
    {
        return getMultiSelectFromListDialog(parent, centerComponent, title, text, objects, objects);
    }

    static public Object[] getMultiSelectFromListDialog(Frame parent, Component centerComponent, String title, String text, Object[] names, Object[] objects)
    {
        StsListSelectionDialog selectDialog = new StsListSelectionDialog(parent, title, text, names);
        selectDialog.setModal(true);
        StsToolkit.centerComponentOnFrame(selectDialog, centerComponent);
        selectDialog.setMultipleSelect();
        selectDialog.setVisible(true);
        int[] selectedIndices = selectDialog.getSelectedIndices();
        if(selectedIndices == null) return null;
        int nSelected = selectedIndices.length;
        if(nSelected == 0) return null;
        Object[] selectedObjects = new Object[nSelected];
        for(int n = 0; n < nSelected; n++)
            selectedObjects[n] = objects[selectedIndices[n]];
        return selectedObjects;
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
            listsPanel.add(label1,
                new GridBagConstraints(0, 0, 1, 1, 0.5, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,
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
        listsPanel.add(list,
            new GridBagConstraints(0, 1, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH,
            new Insets(0, 0, 0, 0), 0, 0));
            Border border = BorderFactory.createEtchedBorder();
        if (button1 != null) button1.addActionListener(
                new java.awt.event.ActionListener()
                {
                    public void actionPerformed(ActionEvent e)
                    {
                        button1_actionPerformed(e);
                    }
                });
        okayButton.addActionListener(
            new java.awt.event.ActionListener()
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
        listsPanel.setBorder(border);
        getContentPane().add(mainPanel);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        setVisible(false);
    }

    protected void cancelButton_actionPerformed(ActionEvent e)
    {
        items = null;
        list.clearSelection();
        setVisible(false);
    }

    protected void button1_actionPerformed(ActionEvent e)
    {
    }

    public int[] getSelectedIndices()
    {
        return list.getSelectedIndices();
    }

    public Object[] getSelectedItems()
    {
        int[] selectedIndices = getSelectedIndices();
        int nSelected = selectedIndices.length;
        Object[] objects = new Object[nSelected];
        for(int n = 0; n < nSelected; n++)
            objects[n] = items[selectedIndices[n]];
        return objects;
    }

    public int getSelectedIndex()
    {
        return list.getSelectedIndex();
    }

    public void setItems(Object[] names)
    {
        items = names;
        list.setListData(names);
    }

    public Object[] getItems() { return items; }

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

    public void setSingleSelect() {list.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION); }
    public void setMultipleSelect() {list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION); }

    public static void main(String[] args)
    {
        String[] names1 = { "Red", "Green", "Blue", "Yellow", "Purple", "Black" };
        StsListSelectionDialog d = new StsListSelectionDialog(new JFrame(), "StsListDialog test", "Double List Example", names1, 1);
        d.setModal(true);
        d.setMultipleSelect();
        d.setVisible(true);
        System.out.print("Selected:");

        Object[] itemsSelected = d.getSelectedItems();
        System.out.print("Selected:");
        for(Object item : itemsSelected)
            System.out.print(" " + item.toString());
        /*
        int[] selectedIndices = d.getSelectedIndices();

        for(int selectedIndex : selectedIndices)
            System.out.print(" " + names1[selectedIndex]);
        */
        
        System.out.println();
        System.exit(0);
    }
}
