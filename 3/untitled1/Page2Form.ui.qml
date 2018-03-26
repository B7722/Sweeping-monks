import QtQuick 2.7
import QtQuick.Controls 2.0

Page {
    width: 600
    height: 400

    title: qsTr("Environment")

    Button {
        id: button
        x: 387
        y: 294
        width: 77
        height: 32
        text: qsTr("Set")
    }

  Label    {
        id: textArea1
        x: 29
        y: 50
        width: 78
        height: 20
        text: qsTr("Lightness")
  }

   Label{
        id: textArea2
        x: 29
        y: 132
        width: 94
        height: 40
        text: qsTr("Position")
   }

    Button {
        id: button1
        x: 476
        y: 294
        width: 63
        height: 32
        text: qsTr("Cancel")
    }

    Label {
        id: textArea3
        x: 29
        y: 219
        text: qsTr("3D")
    }

    TextField {
        id: textField
        x: 129
        y: 50
        width: 110
        height: 32
    }

    Label {
        id: textArea4
        x: 245
        y: 56
        width: 78
        height: 20
        text: qsTr("%")
    }

    ComboBox {
        id: comboBox
        x: 129
        y: 136
        width: 168
        height: 32
    }

    CheckBox {
        id: checkBox2
        x: 129
        y: 219
        text: qsTr("Yes")
    }

    CheckBox {
        id: checkBox3
        x: 237
        y: 219
        text: qsTr("No")
    }

    Image {
        id: image
        x: 0
        y: 0
        width: 593
        height: 350
        z: -1
        source: "images/2.jpg"
    }
}
