package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsTypeDefinePanel extends JPanel implements ActionListener, ListSelectionListener
{
    private StsLithTypesWizard wizard;
    private StsTypeDefine wizardStep;
    private StsModel model = null;

    private String typeName = null;
    private Color typeColor = null;
    private byte[] typeTexture = null;

    JPanel jPanel1 = new JPanel();
    JList typeList = new JList();
    JLabel jLabel1 = new JLabel();
    TitledBorder titledBorder1;
    JPanel jPanel2 = new JPanel();
    JLabel jLabel2 = new JLabel();
    JButton newTypeBtn = new JButton();
    JLabel jLabel3 = new JLabel();
    JLabel jLabel4 = new JLabel();
    JLabel jLabel5 = new JLabel();
    JTextField nameText = new JTextField();
    JTextField textureText = new JTextField();
    JButton textureBtn = new JButton();
    JButton acceptBtn = new JButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    StsColorComboBoxFieldBean colorCombo = new StsColorComboBoxFieldBean();
    StsColor[] colors = null;

    public StsTypeDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsLithTypesWizard)wizard;
        this.wizardStep = (StsTypeDefine)wizardStep;
        this.model = wizard.getModel();

        try
        {
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        colors = model.getSpectrum("Basic").getStsColors();
        colorCombo.initializeColors(this, "typeColor", "Color:", colors);
        colorCombo.setValueObject(Color.red);

        StsObjectRefList types = wizard.getSelectedLibrary().getTypes();
        typeList.setListData(types.getArrayList());
    }

    void jbInit() throws Exception
    {
    titledBorder1 = new TitledBorder("");
    this.setLayout(gridBagLayout3);
    jPanel1.setBorder(BorderFactory.createEtchedBorder());
    jPanel1.setLayout(gridBagLayout2);
    jLabel1.setFont(new java.awt.Font("Dialog", 3, 11));
    jLabel1.setHorizontalAlignment(SwingConstants.LEFT);
    jLabel1.setText("Types");
    typeList.setBorder(BorderFactory.createLoweredBevelBorder());
    jPanel2.setBorder(BorderFactory.createEtchedBorder());
    jPanel2.setDebugGraphicsOptions(0);
    jPanel2.setLayout(gridBagLayout1);
    jLabel2.setText("Selected Type Details");
    jLabel2.setHorizontalAlignment(SwingConstants.LEFT);
    jLabel2.setFont(new java.awt.Font("Dialog", 3, 11));
    newTypeBtn.setText("New");
    jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
    jLabel3.setText("Name:");
    jLabel5.setText("Texture:");
    jLabel5.setHorizontalAlignment(SwingConstants.RIGHT);
    nameText.setText("type name");
    textureText.setText("type texture");
    textureBtn.setText("...");
    acceptBtn.setText("Accept");
    acceptBtn.setEnabled(false);

    typeList.addListSelectionListener(this);
    newTypeBtn.addActionListener(this);
    textureBtn.addActionListener(this);
    acceptBtn.addActionListener(this);

    textureBtn.setEnabled(false);
    textureText.setEnabled(false);
    nameText.setEnabled(false);
    colorCombo.setEditable(false);

    this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 6, 3), 0, -2));
    jPanel1.add(typeList,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 0, 0), 132, 101));
    jPanel1.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 4, 0, 0), 100, 4));
    jPanel1.add(jPanel2,  new GridBagConstraints(1, 1, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 9, 9, 5), 6, -1));
    jPanel2.add(jLabel3,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
    jPanel2.add(nameText,  new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0));
    jPanel1.add(jLabel2,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 9, 0, 5), 122, 8));
    jPanel1.add(newTypeBtn,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 35, 9, 32), 0, 0));
    jPanel2.add(colorCombo,  new GridBagConstraints(0, 1, 3, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 10, 5, 0), 0, 0));
    jPanel2.add(jLabel5,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
    jPanel2.add(textureText,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0));
    jPanel2.add(textureBtn,  new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
    jPanel2.add(acceptBtn,  new GridBagConstraints(1, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        if(source == textureBtn)
        {
            ;
        }
        else if(source == newTypeBtn)
        {
            nameText.setText("");
            colorCombo.setSelectedIndex(0);
            textureText.setText("");
            acceptBtn.setEnabled(true);
            nameText.setEnabled(true);
            colorCombo.setEditable(true);
        }
        else if(source == textureBtn)
        {
            ;
        }
        else if(source == acceptBtn)
        {
            acceptBtn.setEnabled(false);
            nameText.setEnabled(false);
            colorCombo.setEditable(false);

            typeList.removeListSelectionListener(this);
            if(nameText.getText() != null)
            {
                wizard.getSelectedLibrary().getType(nameText.getText(),new StsColor(colorCombo.getStsColor()));
                typeList.removeAll();
                StsObjectRefList types = wizard.getSelectedLibrary().getTypes();
                typeList.setListData(types.getArrayList());
                typeList.setSelectedIndex(types.getSize()-1);
            }
            typeList.addListSelectionListener(this);
        }

    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        if(source == typeList)
        {
            StsObjectRefList types = wizard.getSelectedLibrary().getTypes();
            StsType currentType = (StsType)types.getElement(typeList.getSelectedIndex());
            nameText.setText(currentType.getName());
            colorCombo.setValueObject(currentType.getStsColor());
//            textureText.setText(currentType.);
        }
    }

    public String getTypeName()
    {
        return typeName;
    }

    public Color getTypeColor()
    {
        return typeColor;
    }
    public void setColor(Color color) { this.typeColor = color; }

    public byte[] getTypeTexture()
    {
        return typeTexture;
    }

}
