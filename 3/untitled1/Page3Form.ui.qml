import QtQuick 2.7
import QtQuick.Controls 2.0
Page {
    width: 600
    height: 400

    title: qsTr("Creator")

    BorderImage {
        id: borderImage
        x: 0
        y: 0
        width: 600
        height: 400
        z: -1
        source: "images/2.jpg"
    }

    TextArea {
        id: textArea
        x: 193
        y: 0
        text: qsTr("Draw Yourself")
        font.pointSize: 19
        horizontalAlignment: Text.AlignHCenter

    }

    MouseArea {
        id: mouseArea
        x: 44
        y: 78
        width: 509
        height: 286
        visible: true
    }
}
