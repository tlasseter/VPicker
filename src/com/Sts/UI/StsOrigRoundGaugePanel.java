
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import javax.swing.*;
import java.awt.*;

public class StsOrigRoundGaugePanel extends JPanel
{

   int targetValue = 0;
   int max_value    = 100; // default to 100 for now
   int value       = 22;
   Thread repaintThread;
   int degreesPerSecond = 1;
   Dimension Size;
   double gaugeWidth;
   double gaugeHeight;
   int    centerX = (int)(gaugeWidth/2.0);
   int    centerY = (int)(gaugeHeight/2.0);
   double zeroAngle = 225.0;
   double maxAngle  = -45;
   double range = zeroAngle - maxAngle;

   public StsOrigRoundGaugePanel(Dimension size)
   {
      Size = size;
      gaugeWidth    = Size.width  * 0.9;
      gaugeHeight = Size.height * 0.9;
      centerX = (int)(gaugeWidth/2.0);
      centerY = (int)(gaugeHeight/2.0);
      setSize(Size);
      setMaximumSize(Size);
      setPreferredSize(Size);
   }

   public void updateValue( int i ){ targetValue = i; }
   public void setMaxValue( int i) { max_value = i; }

   private void paintTheBackground( Graphics g)
   {
      g.setColor(this.getBackground());
      g.fillRect(0, 0, Size.width, Size.height);
      g.setColor(Color.darkGray);
      g.fillOval(0, 0, (int)gaugeWidth, (int)gaugeHeight);
      g.setColor(Color.lightGray);
      g.fillOval((int)(gaugeWidth*.01f), (int)(gaugeHeight*.01f), (int)(gaugeWidth-(gaugeWidth*.02f)), (int)(gaugeHeight-gaugeHeight*.02f));
      // now the lines and the arcs on the gauge
      g.setColor(Color.lightGray);
      g.drawLine(centerX, centerY, (int) gaugeWidth -25, (int)gaugeHeight-25);
      g.drawLine(centerX, centerY, 23, (int)gaugeHeight-25);
      g.setColor(Color.darkGray);
      g.fillArc(10, 10, (int)gaugeWidth-20, (int)gaugeHeight-20, -45, 270);
      g.setColor(Color.blue);
      g.drawArc( 10, 10, (int)gaugeWidth-20, (int)gaugeHeight-20, -45, 270);
      // this red line doesn't add much to the gauge so....
      //g.setColor( Color.red);
      //g.drawArc( 10, 10, (int)gaugeWidth-20, (int)gaugeHeight-20, -45, 45 );
      g.setColor(Color.white);
      g.drawString("0%", centerX - 45, centerY + 55);
      g.drawString("50%", centerX - 5, centerY - 50);
      g.drawString("100%", centerX + 15, centerY + 55);

   }

   public void paintComponent(Graphics g)
   {
      Color oldColor;
      oldColor = g.getColor();
      paintTheBackground( g);

      g.setColor(Color.red);
      int x1 = centerX,
          x2 = x1,
          y1 = centerY,
          y2 = y1;
      double angleToUse = zeroAngle - 1.0 * range *( value * 1.0 / max_value * 1.0);
      x2 += (int)( Math.cos(Math.toRadians(angleToUse))*centerX);
      y2 -= (int)( Math.sin(Math.toRadians(angleToUse))*centerY);
      g.drawLine(x1, y1, x2, y2 );
      //g.setColor(Color.white);
      g.drawString(""+ value, centerX - 10, centerY + 30);
      g.setColor(oldColor);
   }

}