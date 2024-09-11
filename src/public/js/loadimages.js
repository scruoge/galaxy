function Cache()
{
  var _imgs = {};
  var _addedImageCount = 0;
  var _loadedImgsCount = 0;
  this.addSpriteSource = function(src) {

    var img = new Image();
    img.onload = function() {
      _loadedImgsCount++;
    };

    img.src = src;
    _imgs[src] = img;
    _addedImageCount++;
  };

  this.getLoadedImagePc = function() {
    return _loadedImgsCount * 100 / _addedImageCount;
  };

  this.getImage = function(src) {
    return _imgs[src];
  };
}

MyCache = new Cache();

MyCache.addSpriteSource("img/stars/small.png");
MyCache.addSpriteSource("img/stars/star_bg2.jpg");

MyCache.addSpriteSource("img/hud.png");
MyCache.addSpriteSource("img/btn_right.png");
MyCache.addSpriteSource("img/btn_down.png");
MyCache.addSpriteSource("img/btn_left.png");
MyCache.addSpriteSource("img/btn_up.png");
MyCache.addSpriteSource("img/galaxy_bg1.png");

function waitImagesLoading()
{
  var pc = MyCache.getLoadedImagePc();
  if (pc < 100)
    setTimeout(waitImagesLoading, 10);
}
