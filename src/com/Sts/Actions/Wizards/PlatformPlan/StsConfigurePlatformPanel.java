package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>e
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsConfigurePlatformPanel extends StsFieldBeanPanel implements ActionListener
{
    private StsPlatform platform;

    StsGroupBox platformPropertiesBox = new StsGroupBox("Platform Description");
    StsStringFieldBean platformNameBean = new StsStringFieldBean();
    StsComboBoxFieldBean platformTypeBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean slotTypeBean = new StsComboBoxFieldBean();

    StsGroupBox slotPropertiesBox = new StsGroupBox("Slot Configuration");
    StsIntFieldBean nRowsBean = new StsIntFieldBean();
    StsFloatFieldBean rowSpacingBean = new StsFloatFieldBean();
    StsIntFieldBean nColsBean = new StsIntFieldBean();
    StsFloatFieldBean colSpacingBean = new StsFloatFieldBean();
    StsDoubleFieldBean xOriginBean = new StsDoubleFieldBean();
    StsDoubleFieldBean yOriginBean = new StsDoubleFieldBean();
    StsFloatFieldBean zKbBean = new StsFloatFieldBean();
    StsFloatFieldBean rotationAngleBean = new StsFloatFieldBean();
    JCheckBox labelsOnChk = new JCheckBox("Show Slot Labels");

    public boolean labelSlots = true;

    StsGroupBox slotGraphicBox = new StsGroupBox("Slot Graphic");
    StsPlatformConfigPanel slotGraphicPanel = new StsPlatformConfigPanel();

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsConfigurePlatformPanel(StsWizard wizard)
    {
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize(StsPlatform platform)
    {
        this.platform = platform;
        setPanelObject(platform);
        slotGraphicPanel.setPlatform(platform);
        platform.setCanvas(slotGraphicPanel);
        slotGraphicPanel.setLabelSlots(labelSlots);
        platform.drawConfiguration();
        slotGraphicPanel.setEnabled(false);
    }

    void jbInit() throws Exception
    {
        setLayout(gridBagLayout1);

        platformNameBean.classInitialize(StsPlatform.class, "name", true, "Name:");
        platformTypeBean.classInitialize(StsPlatform.class, "platformTypeByString", "Platform Type:", StsPlatform.PLATFORM_TYPE_STRINGS);
        slotTypeBean.classInitialize(StsPlatform.class, "slotLayoutByString", "Slot Configuration:", StsPlatform.SLOT_LAYOUT_STRINGS);

        nColsBean.classInitialize(StsPlatform.class, "nCols", true, "Slot Columns:");
        colSpacingBean.classInitialize(StsPlatform.class, "colSpacing", true, "Column Spacing:");
        nRowsBean.classInitialize(StsPlatform.class, "nRows", true, "Slot Rows:");
        rowSpacingBean.classInitialize(StsPlatform.class, "rowSpacing", true, "Row Spacing:");
        xOriginBean.classInitialize(StsPlatform.class, "xOrigin", true, "X Origin:");
        yOriginBean.classInitialize(StsPlatform.class, "yOrigin", true, "Y Origin:");
        zKbBean.classInitialize(StsPlatform.class, "zKB", true, "KB height:");
        rotationAngleBean.classInitialize(StsPlatform.class, "rotationAngle", true, "Rotation Angle:");
        labelsOnChk.addActionListener(this);
        labelsOnChk.setSelected(true);
//        reLabelSlotsBtn.addActionListener(this);

        platformPropertiesBox.add(platformNameBean);
        platformPropertiesBox.add(platformTypeBean);
        platformPropertiesBox.add(slotTypeBean);
        add(platformPropertiesBox);

        slotPropertiesBox.add(nColsBean);
        slotPropertiesBox.add(colSpacingBean);
        slotPropertiesBox.add(nRowsBean);
        slotPropertiesBox.add(rowSpacingBean);
        slotPropertiesBox.add(xOriginBean);
        slotPropertiesBox.add(yOriginBean);
        slotPropertiesBox.add(zKbBean);
        slotPropertiesBox.add(rotationAngleBean);
        slotPropertiesBox.add(labelsOnChk);
//        slotBeanPanel.add(reLabelSlotsBtn);
        add(slotPropertiesBox);

        slotGraphicPanel.setBackground(Color.WHITE);
        slotGraphicPanel.setSize(300,200);
        slotGraphicBox.add(slotGraphicPanel, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 300, 200));
        add(slotGraphicBox);
    }

    public void setXOrigin(double x)
    {
        xOriginBean.setValue(x);
    }

    public void setYOrigin(double y)
    {
        yOriginBean.setValue(y);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == labelsOnChk)
        {
            slotGraphicPanel.setLabelSlots(labelsOnChk.isSelected());
            platform.drawConfiguration();
        }
//        else if(source == reLabelSlotsBtn)
//        {
//            new StsMessage(wizard.frame, StsMessage.WARNING, "Function not currently available");
//        }
    }

    public void setLabelSlots(boolean val)
    {
        slotGraphicPanel.setLabelSlots(val);
        platform.drawConfiguration();
    }
    public boolean getLabelSlots() { return slotGraphicPanel.getLabelSlots(); }

    public void setEditable(boolean enabled)
    {
        platformPropertiesBox.setEnabled(false);
        slotPropertiesBox.setEnabled(false);
    }

/*
    static class Demo extends JPanel
    {
        private BufferedImage bimg;
        private StsPlatform platform;
        private Ellipse2D ellipse = null;
        private Line2D line = null;
        private boolean labelSlots = false;

        public Demo()
        {
            setBackground(Color.white);
        }

        public void setPlatform(StsPlatform platform)
        {
            this.platform = platform;
        }

        public void setLabelSlots(boolean val)
        {
            labelSlots = val;
        }
        public boolean getLabelSlots() { return labelSlots; }

        public void drawDemo(int w, int h, Graphics2D g2)
        {
            if(platform == null)
                return;
            g2.setBeachballColors(Color.blue);
            Font font = g2.getFont();
            FontRenderContext frc = g2.getFontRenderContext();
            TextLayout tl = new TextLayout(platform.getName(), font, frc);
            tl.draw(g2, (float) (w/2-tl.getBounds().getWidth()/2),(float) (tl.getAscent()+tl.getDescent()));

            if(platform.checkValidity())
            {
                if(platform.getSlotLayoutByString().equals(platform.SLOT_LAYOUT_STRINGS[platform.RECTANGULAR]))
                {
                    float dw = w * 0.6f;
                    float dh = h * 0.6f;
                    float dwi = dw / platform.getNCols();
                    float dhi = dh / platform.getNRows();
                    if(dwi > dhi)
                        dwi = dhi;
                    else
                        dhi = dwi;
                    if(platform.getColSpacing() > platform.getRowSpacing())
                        dhi = dhi * (platform.getRowSpacing() / platform.getColSpacing());
                    else
                        dwi = dwi * (platform.getColSpacing() / platform.getRowSpacing());
                    float sw = w * 0.2f;
                    float sh = h * 0.15f;

                    // Origin Text
                    g2.setBeachballColors(Color.blue);
                    String originStg = new String(String.valueOf("(" + platform.getXOrigin() + ", " + String.valueOf(platform.getYOrigin()) + ")"));
                    tl = new TextLayout(originStg, font, frc);
                    tl.draw(g2, (float)(w/2-(tl.getBounds().getWidth()/2.0f)), h + (tl.getAscent()+tl.getDescent()*2));

                    // Vertical North line
                    g2.setBeachballColors(Color.red);
                    line = new Line2D.Float(w *0.05f, h*0.10f, w *0.05f, h*0.05f );
                    g2.draw(line);
                    line = new Line2D.Float(w *0.06f, h*0.10f, w *0.06f, h*0.05f );
                    g2.draw(line);
                    line = new Line2D.Float(w *0.045f, h*0.07f, w *0.055f, h*0.03f );
                    g2.draw(line);
                    line = new Line2D.Float(w *0.065f, h*0.07f, w *0.055f, h*0.03f );
                    g2.draw(line);
                    tl = new TextLayout("North", font, frc);
                    tl.draw(g2, (float)(w*0.055f -tl.getBounds().getWidth()/2.0f) , (float)(h*.10f + (tl.getAscent()+tl.getDescent()*2)));

                    // Create Rotation Transform
                    float centerX = sw + dw/2.0f;
                    float centerY = h - sh - dh/2;
                    AffineTransform aT = g2.getTransform();
                    AffineTransform at = AffineTransform.getTranslateInstance(centerX, centerY);
                    at.rotate(Math.toRadians(platform.getRotationAngle()));
                    g2.transform(at);

                    // Draw center of pad
                    g2.setBeachballColors(Color.blue);
                    line = new Line2D.Float(2, 0, 5, 0);
                    g2.draw(line);
                    line = new Line2D.Float(-2, 0, -5, 0);
                    g2.draw(line);
                    line = new Line2D.Float(0, 2, 0, 5 );
                    g2.draw(line);
                    line = new Line2D.Float(0, -2, 0, -5 );
                    g2.draw(line);

                    // Rotating angle line
                    line = new Line2D.Float(0, -dh/2, 0, -(dh/2 + h*0.075f));
                    g2.draw(line);
                    line = new Line2D.Float(w*0.01f, -dh/2, w*0.01f, -(dh/2 + h*0.075f));
                    g2.draw(line);
                    line = new Line2D.Float(-w*0.01f, -(dh/2+h*0.055f), w*0.005f, -(dh/2 + h*0.10f));
                    g2.draw(line);
                    line = new Line2D.Float(w*0.02f, -(dh/2+h*0.055f), w*0.005f, -(dh/2 + h*0.10f));
                    g2.draw(line);

                    // Rotation Angle Text
                    g2.setBeachballColors(Color.blue);
                    tl = new TextLayout(String.valueOf(platform.getRotationAngle()), font, frc);
                    tl.draw(g2, (float) -tl.getBounds().getWidth()/2.0f, -(dh/2 + h*0.10f)-tl.getDescent());

                    // Slots
                    g2.setBeachballColors(Color.black);
                    float startX = -4f;
                    if(platform.getNCols() > 1)
                        startX = -(dwi*(platform.getNCols()-1))/2.0f - 4f;
                    float startY = -4f;
                    if(platform.getNRows() > 1)
                        startY = (dhi*(platform.getNRows()-1))/2.0f - 4f;

                    for (int j = 0; j < platform.getNRows(); j++)
                        for (int i = 0; i < platform.getNCols(); i++)
                        {
                            ellipse = new Ellipse2D.Float(startX + (dwi * i), startY + (-dhi * j), 8, 8);
                            g2.draw(ellipse);
                            if(labelSlots)
                            {
                                tl = new TextLayout(String.valueOf((j*platform.getNCols()) + i + 1), font, frc);
                                tl.draw(g2, (float) startX + (dwi * i), startY + (-dhi * j));
                            }
                        }
                    g2.transform(aT);
                }
                else
                {
                    // Circular Layout
                    g2.setBeachballColors(Color.red);
                    tl = new TextLayout("Circular Layout is not currently supportted", font, frc);
                    tl.draw(g2, (float)(w/2-tl.getBounds().getWidth()/2),(float)(h/2));
                }
            }
            else
            {
                TextLayout t2 = new TextLayout("Incomplete Definition", font, frc);
                g2.setBeachballColors(Color.red);
                t2.draw(g2, (float) (w / 2 - t2.getBounds().getWidth() / 2), (float) (h /2));
            }
        }

        public Graphics2D createGraphics2D(int w, int h)
        {
            Graphics2D g2 = null;
            if (bimg == null || bimg.getWidth() != w || bimg.getHeight() != h)
                bimg = (BufferedImage) createImage(w, h);

            g2 = bimg.createGraphics();
            g2.setBackground(getBackground());
            g2.clearRect(0, 0, w, h);
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_RENDERING,
                                RenderingHints.VALUE_RENDER_QUALITY);
            return g2;
        }

        public void paint(Graphics g)
        {
            Dimension d = getSize();
            Graphics2D g2 = createGraphics2D(d.width, d.height);
            drawDemo(d.width, d.height, g2);
            g2.dispose();
            g.drawImage(bimg, 0, 0, this);
        }
    }
    */
}
