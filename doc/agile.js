/* $Id$ */

/* Javascript collection for AGILE */

/* Function to create the flexible table of contents */
function createFlexiToC()
{
  var headerNodes = $(':header');

  // no ToC for less than two elements
  if (headerNodes.length < 2)
    return false;

  // --- now it is time to create the ToC ---
  var flexiToC = $('<div id="flexitoc"></div>');

  // add caption with toggle event
  var tocCaption = $('<div class="toc-caption">show contents</div>');
  tocCaption.toggle(function() {
      $('#flexitoc div.toc-caption').html('hide contents');
      $('#flexitoc div.toc-listing').slideDown('slow');
  }, function() {
      $('#flexitoc div.toc-caption').html('show contents');
      $('#flexitoc div.toc-listing').slideUp('slow');
  });
  flexiToC.append(tocCaption);

  // add content listing
  var tocListing = $('<div class="toc-listing" style="display:none;"></div>');
  headerNodes.each(function(idx) {
    var $this = $(this);
    var anchor = $this.children('a:first-child');
    var text = anchor.length > 0 ? anchor.text() : $this.text();
    /* create the link -> the heading must get a unique ID first, if it does
       not already have one */
    var cid = $this.attr('id') || ('toc_link' + idx);
    $this.attr('id', cid);
    var class = 'page';
    if (this.tagName == 'H1') {
      cid = 'top';
      text = 'Top';
    } else {
      var eCnt = parseInt(this.tagName.substr(1));
      if (eCnt >= 3) {
	class += ' ';
	for (var i=0; i<(eCnt-3); ++i)
	  class += 'e';
	class += 'indent';
      }
    }
    tocListing.append('<a class="'+class+'" href="#'+cid+'">'+text+'</a>');
  });
  flexiToC.append(tocListing);

  // add toc to dom tree
  $('#TOC').html(flexiToC);

  // set scroll callback, and call for init
  $(window).scroll(function() {
    var $toc = $('#TOC');
    var stdOffset = $toc.data('initialOffset');
    if (typeof(stdOffset) === 'undefined') {
      // initialization
      $toc.data('initialOffset', $toc.offset());
      return;
    }
    var newOffsetTop = $(window).scrollTop() + 5;
    newOffsetTop = stdOffset.top > newOffsetTop ? stdOffset.top : newOffsetTop;
    $toc.css('top', parseInt(newOffsetTop));
  }).scroll();
  return true;
}


// jQuery onLoad hook
$(createFlexiToC);
$(function() {
  // init unittest links
  if (typeof(UNITTEST_EXISTS) !== 'undefined' && UNITTEST_EXISTS) {
    $('.unittestlink .unittestexists').show();
    $('.unittestlink .unittestmissing').hide();
  } else {
    $('.unittestlink .unittestexists').hide();
    $('.unittestlink .unittestmissing').show();
  }
});

/* End of $Id$ */
