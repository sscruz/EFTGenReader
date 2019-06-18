import os

from HTMLGenerator import *
# Create an index.html file for a web-directory with .png files

FCN_STR = """
function myFunction() {
    // Declare variables
    var input, filter, u, li, a, i, txtValue;
    input = document.getElementById('myInput');
    filter = input.value.toUpperCase();

}
"""

STYLE_STR = """
.image {
    float:left; margin: 5px; clear:justify;
    font-size: 10px; font-family: Verdana, Arial, sans-serif;
    text-align: center;
}

.topnav {
    overflow: hidden;
    background-color: #e9e9e9;
}

.topnav a {
    float: left;
    display: block;
    color: black;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none;
}

.topnav a:hover {
    background-color: #ddd;
    color: black;
}

.topnav a.active {
    background-color: #2196F3;
    color: white;
}

.topnav input[type=text] {
    float: right;
    padding: 6px;
    margin-top: 8px;
    margin-right: 16px;
    border: none;
    font-size: 17px;
}
"""

def getImages(tar_dir,file_type='png'):
    # Note: Will break if a filename has more then 1 . in the name
    fnames = []
    for out in sorted(os.listdir(tar_dir)):
        fpath = os.path.join(tar_dir,out)
        if (os.path.isdir(fpath)):
            continue
        f,ftype = out.split(".")
        if ftype != file_type:
            continue
        fnames.append(out)
    return fnames

def make_html(tar_dir):
    home_dir = os.getcwd()
    if not os.path.exists(tar_dir):
        print "Target directory does not exists: %s" % (tar_dir)
        return

    os.chdir(tar_dir)

    my_html = HTMLGenerator()

    meta_tag = MetaTag(); my_html.addHeadTag(meta_tag)
    meta_tag.addAttributes(charset='UTF-8')

    style_tag = StyleTag(); my_html.addHeadTag(style_tag)
    style_tag.setContent(STYLE_STR)
    style_tag.addAttributes(type='text/css')

    #topnav = DivisionTag()
    #topnav.addAttributes(cls='topnav')
    #input_tag = InputTag()
    #input_tag.addAttributes(type='text',placeholder='Search...',id='myInput')
    #my_html.addBodyTag(topnav)
    #topnav.addTag(input_tag)

    image_files = getImages(tar_dir)
    for fname in image_files:
        image_name,ftype = fname.split(".")

        div_tag   = DivisionTag(); my_html.addBodyTag(div_tag)
        image_tag = ImgTag()
        text_div  = DivisionTag()
        link_tag  = HyperLinkTag(link_location="./%s" % (fname),link_name='')

        # This ensures the pretty_print setting gets inherited properly
        div_tag.addTag(link_tag); div_tag.addTag(text_div)
        link_tag.addTag(image_tag)


        image_tag.addAttributes(width=355,height=229,border=0,src="./%s" % (fname))

        link_tag.addAttributes(target='_blank')
        
        text_div.addAttributes(style='width:355px',id='imgName')
        text_div.setContent(image_name)

        #div_tag.addAttribute('class','image')
        div_tag.addAttributes(cls='image')


    print my_html.dumpHTML()
    my_html.saveHTML(f_name='index.html',f_dir=tar_dir)


    os.chdir(home_dir)

def main():
    web_area = "/afs/crc.nd.edu/user/a/awightma/www/"
    if len(sys.argv) == 2:
        dir_name = os.path.join("eft_stuff/misc/gen_plots",sys.argv[1])
    else:
        dir_name = "eft_stuff/tmp"
    fpath = os.path.join(web_area,dir_name)
    make_html(fpath)

if __name__ == "__main__":
    main()