package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWeighPointDialog extends JDialog implements ActionListener, ChangeListener
{
    private ButtonGroup btnGroup1 = new ButtonGroup();
    public StsModel model = null;
    public JPanel jPanel1 = new JPanel();
    public JPanel buttonPanel = new JPanel();
    public JButton newBtn = new JButton("New");
    public JButton exitBtn = new JButton("Exit");
    Font defaultFont = new Font("Dialog",0,11);
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JTextArea notesArea = new JTextArea();
    StsWeighPoint wp = null;
    DateFormat dateFormat = DateFormat.getDateTimeInstance();

    JPanel jPanel2 = new JPanel();
    JTextArea newNote = new JTextArea();
    TitledBorder titledBorder1;
    JSpinner noteSpinner = new JSpinner();
    JTextField nameTxt = new JTextField();
    JLabel jLabel1 = new JLabel();
    SpinnerListModel spinnerModel = null;
    String[] list = null;
  GridBagLayout gridBagLayout4 = new GridBagLayout();
  GridBagLayout gridBagLayout5 = new GridBagLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
  JButton showAllBtn = new JButton();
  GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsWeighPointDialog(StsModel model, StsWeighPoint weighPoint, boolean modal)
    {
        super(model.win3d ,"WayPoint Viewer / Editor", modal);
        this.setLocationRelativeTo(model.win3d);
        this.model = model;
        wp = weighPoint;
        initialize();
        this.setSize(new Dimension(250,400));

        try
        {
            jbInit();
            pack();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        titledBorder1 = new TitledBorder("");
        this.setModal(false);
        this.getContentPane().setLayout(gridBagLayout5);

        newBtn.setFont(defaultFont);
        newBtn.setToolTipText("Add a message");
        newBtn.setActionCommand("New");
        newBtn.setText("Add");
        newBtn.addActionListener(this);
        exitBtn.setFont(defaultFont);
        exitBtn.addActionListener(this);

        buttonPanel.setBorder(BorderFactory.createEtchedBorder());
        buttonPanel.setLayout(gridBagLayout3);

        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setToolTipText("");
        jPanel2.setLayout(gridBagLayout2);
        newNote.setBorder(BorderFactory.createLoweredBevelBorder());
        newNote.setMinimumSize(new Dimension(325, 65));
        newNote.setPreferredSize(new Dimension(325, 65));
        newNote.setToolTipText("Enter your note here");
        newNote.setLineWrap(true);
        newNote.setWrapStyleWord(true);
        notesArea.setBackground(Color.lightGray);
        notesArea.setFont(new java.awt.Font("DialogInput", 0, 11));
        notesArea.setBorder(BorderFactory.createLineBorder(Color.black));
        notesArea.setCaretPosition(0);
        notesArea.setEditable(false);
        notesArea.setColumns(0);
        notesArea.setLineWrap(true);
        notesArea.setWrapStyleWord(true);

        nameTxt.setBackground(Color.lightGray);
        nameTxt.setBorder(BorderFactory.createEtchedBorder());
        nameTxt.setMinimumSize(new Dimension(75, 19));
        nameTxt.setToolTipText("");
        nameTxt.setMargin(new Insets(1, 1, 1, 1));
        nameTxt.setEditable(true);
        nameTxt.setText("WayPointName");
        nameTxt.setHorizontalAlignment(SwingConstants.LEADING);
        jLabel1.setHorizontalAlignment(SwingConstants.LEFT);
        jLabel1.setText("New Note");
        noteSpinner.setFont(new java.awt.Font("Dialog", 1, 10));
        noteSpinner.setOpaque(false);
        noteSpinner.setToolTipText("Display all notes after a certain time");
        noteSpinner.setFont(new java.awt.Font("DialogInput", 0, 11));
        jPanel1.setLayout(gridBagLayout4);
        showAllBtn.setToolTipText("Show all the message");
        showAllBtn.setText("Show All");
        showAllBtn.addActionListener(this);
        buttonPanel.add(exitBtn,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 61, 3, 5), 0, 0));
        buttonPanel.add(newBtn,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 3, 0), 0, 0));
        buttonPanel.add(showAllBtn,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 64, 3, 0), 14, 1));
        jPanel1.add(jPanel2,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 0, 0), -159, 89));
        jPanel2.add(newNote,   new GridBagConstraints(0, 3, 2, 1, 1.0, 0.5
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 3, 5, 4), 2, 28));
        jPanel2.add(noteSpinner,     new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 4), 171, 0));
        jPanel2.add(nameTxt,      new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 3, 0, 3), 136, 0));
        jPanel2.add(jLabel1,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 3, 0, 0), 79, 2));
        jPanel2.add(jScrollPane1,   new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 3, 0, 4), 0, 0));
        jScrollPane1.getViewport().add(notesArea, null);
        jPanel1.add(buttonPanel,     new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 6, 0), 0, 0));

        this.getContentPane().add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 6, 4), 0, 0));
    nameTxt.setEditable(false);

    }

    private void initialize()
    {

        list = null;
        spinnerModel = null;

        // Update Name
        nameTxt.setText(wp.getName());
        newNote.setText("");

        // Add the notes to the text area
        if(wp.getNotes() == null)
        {
            notesArea.setText("None");
            return;
        }
        else
        {
            buildSpinner();
            setCurrentNote(wp.getNotes().getSize()-1);
        }
    }

    private void buildSpinner()
    {
        long date = 0L;
        int i1 = 0;

        list = new String[wp.getNotes().getSize()];
        for(int i=0; i<wp.getNotes().getSize(); i++)
        {
            date = ((StsNote)wp.getNotes().getElement(i)).timeStamp;
            i1 = i+1;
            list[i] = new String("Note #" + i1 + " " + dateFormat.format(new Date(date)));
        }
        spinnerModel = new SpinnerListModel(list);
        noteSpinner.setModel(spinnerModel);
    }

    public void setCurrentNote(int idx)
    {
        long date = 0L;
        String noteString = null;
        StsNote note = null;

        noteSpinner.removeChangeListener(this);
        if(wp.getNotes() == null)
        {
            new StsMessage(model.win3d, StsMessage.WARNING,
                           "No notes exist for this waypoint.");
            return;
        }

        for(int i=idx; i<wp.getNotes().getSize(); i++)
        {
            date = ( (StsNote) wp.getNotes().getElement(i)).timeStamp;
            note = (StsNote) wp.getNotes().getElement(i);
            if(noteString != null)
                noteString = new String(noteString + "\n\nAuthor: " + note.author + "\nDate: " +
                                    dateFormat.format(new Date(date))
                                    + "\n     " + note.note);
            else
                noteString = new String("Author: " + note.author + "\nDate: " +
                                    dateFormat.format(new Date(date))
                                    + "\n     " + note.note);

        }
        notesArea.setText(noteString);
        noteSpinner.setValue(list[idx]);
        noteSpinner.addChangeListener(this);
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        if(source == noteSpinner)
        {
            String item = (String)noteSpinner.getValue();
            for(int i=0; i<list.length; i++)
                if(list[i].equalsIgnoreCase(item))
                    setCurrentNote(i);
        }
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == newBtn)
        {
            // Add a new note
            if(newNote.getText() != null)
            {
                wp.addNote(newNote.getText());
                initialize();
            }
            else
                StsMessageFiles.infoMessage("No note has been entered.");
        }
        else if(source == showAllBtn)
        {
            setCurrentNote(0);
        }
        else if(source == exitBtn)
        {
            // Exit this dialog
            this.dispose();
        }
    }
}
