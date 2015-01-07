package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsArbitraryLineDefinePanel extends StsFieldBeanPanel implements ActionListener
{
    StsArbitraryLineWizard wizard;
    StsArbitraryLineDefine step;

    static final int START = 0;
    static final int CENTER = 1;
    static final int END = 2;

    StsGroupBox defineBox = new StsGroupBox("Define Hinge Point");
    public StsToggleButton createButton;
    public StsToggleButton deleteButton;
    public StsToggleButton moveButton;
    public StsCheckbox wellSnapChk;

    StsGroupBox hingeBox = new StsGroupBox("Hinge Parameters");
    JLabel listLbl = new JLabel("Hinges:");
    JComboBox hingeList = new JComboBox();
    JLabel rowLbl = new JLabel("InLine:");
    JTextField rowTxt = new JTextField();
    JLabel colLbl = new JLabel("CrossLine:");
    JTextField colTxt = new JTextField();
    JLabel startLbl = new JLabel("Start:");
    JTextField startTxt = new JTextField();
    JLabel endLbl = new JLabel("End:");
    JTextField endTxt = new JTextField();

    StsGroupBox rotateBox = new StsGroupBox("Rotation Parameters");
    JRadioButton startRadio = new JRadioButton();
    JRadioButton centerRadio = new JRadioButton();
    JRadioButton endRadio = new JRadioButton();
    ButtonGroup rotateRadioGrp = new ButtonGroup();
    JLabel jLabel1 = new JLabel();
    JTextField rotateAngleTxt = new JTextField();
    JLabel jLabel2 = new JLabel();

    GridBagLayout gridBagLayout1 = new GridBagLayout();

    StsLine hinge = null;
    int row, col, start, end;

    static final String createButtonTip = "Create new hinge point";
    static final String deleteButtonTip = "Select and delete a hinge point";
    static final String moveButtonTip = "Select and move an existing hinge point";
    static final String wellSnapTip = "Snap hinge to the nearest well";
    static final String startRadioTip = "Rotate about the first point defined in line";
    static final String centerRadioTip = "Rotate about the center point of the line";
    static final String endRadioTip = "Rotate about the last point defined in line";
  JButton rotateBtn = new JButton();
  GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsArbitraryLineDefinePanel(StsWizard wizard, StsWizardStep step)
    {
        this.wizard = (StsArbitraryLineWizard)wizard;
        this.step = (StsArbitraryLineDefine)step;

        try
        {
            jbInit();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        setLayout(gridBagLayout1);
        rotateBox.setLayout(gridBagLayout2);
        startRadio.setText("Start");
        centerRadio.setText("Center");
        endRadio.setText("End");
        startRadio.setToolTipText(startRadioTip);
        centerRadio.setToolTipText(centerRadioTip);
        centerRadio.setSelected(true);
        endRadio.setToolTipText(endRadioTip);

        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Rotate About:");
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel2.setText("Angle:");
        rotateBtn.setText("Rotate");
        rotateAngleTxt.setSelectionEnd(0);
    rotateAngleTxt.setText("0.0");
    rotateAngleTxt.setHorizontalAlignment(SwingConstants.RIGHT);
    add(hingeBox,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 174, 41));

        createButton = new StsToggleButton("Create", createButtonTip, this, "createHinge");
        deleteButton = new StsToggleButton("Delete", deleteButtonTip, this, "deleteHinge");
        moveButton = new StsToggleButton("Move", moveButtonTip, this, "moveHinge");

        wellSnapChk = new StsCheckbox();
        wellSnapChk.setText("Snap to Well");
        wellSnapChk.setSelected(false);
        wellSnapChk.setEnabled(false);

        defineBox.add(createButton, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 50, 10));
        defineBox.add(moveButton, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 50, 10));
        defineBox.add(deleteButton, new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 50, 10));
        defineBox.add(wellSnapChk, new GridBagConstraints(0, 1, 1, 3, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
        this.add(rotateBox,   new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 3, 5), 263, 0));
        this.add(defineBox,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(7, 5, 0, 5), 91, 0));

        hingeBox.setLayout(new GridBagLayout());
        rowTxt.setEnabled(false);
        rowTxt.setBorder(BorderFactory.createEtchedBorder());
        rowTxt.setBackground(Color.lightGray);
        colTxt.setEnabled(false);
        colTxt.setBorder(BorderFactory.createEtchedBorder());
        colTxt.setBackground(Color.lightGray);
        startTxt.setEnabled(false);
        startTxt.setBorder(BorderFactory.createEtchedBorder());
        startTxt.setBackground(Color.lightGray);
        endTxt.setEnabled(false);
        endTxt.setBorder(BorderFactory.createEtchedBorder());
        endTxt.setBackground(Color.lightGray);

        hingeBox.add(listLbl, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
           GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
        hingeBox.add(hingeList, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
           GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 100, 0));
        hingeBox.add(rowLbl, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        hingeBox.add(rowTxt, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));
        hingeBox.add(colLbl, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        hingeBox.add(colTxt, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));
        hingeBox.add(startLbl, new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        hingeBox.add(startTxt, new GridBagConstraints(1, 3, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));
        hingeBox.add(endLbl, new GridBagConstraints(0, 4, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        hingeBox.add(endTxt, new GridBagConstraints(1, 4, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));
        rotateRadioGrp.add(endRadio);
        rotateRadioGrp.add(centerRadio);
        rotateRadioGrp.add(startRadio);
    rotateBox.add(endRadio,  new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 10, 0, 23), 14, -6));
    rotateBox.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(7, 29, 0, 0), 20, 3));
    rotateBox.add(startRadio,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 6, 0, 0), 11, -6));
    rotateBox.add(centerRadio,  new GridBagConstraints(2, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 11, 0, 0), 13, -6));
    rotateBox.add(rotateAngleTxt,  new GridBagConstraints(1, 1, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(11, 9, 0, 0), 32, 2));
    rotateBox.add(rotateBtn,  new GridBagConstraints(3, 1, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 14, 0, 32), 28, 0));
    rotateBox.add(jLabel2,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(11, 44, 0, 0), 40, 7));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == rotateBtn)
        {
            activateRotation(new Float(rotateAngleTxt.getText()).floatValue());
        }
    }

    // Rotate the line
    public void activateRotation(float angle)
    {

    }

    public int getRotateOrigin()
    {
        if(startRadio.isSelected())
        {
            return START;
        }
        else if(centerRadio.isSelected())
        {
            return CENTER;
        }
        else
            return END;
    }

    public float getRotateAngle()
    {
        return new Float(rotateAngleTxt.getText()).floatValue();
    }

    public void initialize()
    {
        int nItems = initializeComboBox();
        createButton.setEnabled(true);
        moveButton.setEnabled(nItems > 0);
        deleteButton.setEnabled(nItems > 0);
    }

    // Enable graphical creation of a hinge
    public void createHinge()
    {
        boolean selected = createButton.isSelected();
        wizard.setSelectHinge(selected);
        // Set the inline, crossline, startz and endz fields.
    }

    // The user has selected an existing hinge
    public void hingeSelected(StsLine hinge)
    {
        setHingeSelected(hinge);
    }

    private int initializeComboBox()
    {
        StsObjectRefList hinges = wizard.getHinges();
        if(hinges == null) return 0;
        int nHinges = hinges.getSize();
        if(nHinges == 0)
        {
            ;
        }
        else
        {
            ;
        }
        return nHinges;
    }

    public void addHingeToComboBox(StsLine hinge)
    {
        return;
    }

    public void setHingeSelected(StsLine hinge)
    {
        if(hinge == null) return;
    }

    // Enable graphical deletion of hinge
    public void deleteHinge()
    {


    }

    // Enable graphical movement of hinge
    public void moveHinge()
    {

    }
}
